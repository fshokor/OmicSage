"""
pipeline/modules/qc/cluster.py
OmicSage — Leiden clustering module

Input  : reduced AnnData (output of reduce.py)
         Required: obsm['X_pca'], obsp['connectivities']
Output : (adata_clustered, metrics_dict)

Resolution selection strategy (in priority order)
--------------------------------------------------
1. best_resolution_override   — manual pin, bypasses all scoring
2. n_clusters_expected        — pick resolution whose cluster count is closest
                                 to the biological prior; use stability as tiebreaker
3. stability plateau          — find where adding resolution stops adding clusters;
                                 use silhouette as tiebreaker
4. silhouette argmax          — fallback when nothing else is specified

Provenance is stored in adata.uns['omicsage_cluster'].
"""

from __future__ import annotations

import importlib.metadata
import logging
import warnings
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from sklearn.metrics import silhouette_score

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

DEFAULT_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
_SILHOUETTE_MAX_CELLS = 10_000


def cluster(
    adata: AnnData,
    resolution_range: list[float] = DEFAULT_RESOLUTIONS,
    best_resolution_override: Optional[float] = None,
    n_clusters_expected: Optional[int] = None,
    pca_key: str = "X_pca",
    connectivities_key: str = "connectivities",
    random_state: int = 0,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Run Leiden clustering at multiple resolutions and auto-select the best one.

    Parameters
    ----------
    adata : AnnData
        Reduced AnnData produced by reduce().  Must contain:
          - obsm[pca_key]               PCA embedding
          - obsp[connectivities_key]    kNN connectivity graph
    resolution_range : list of float
        Leiden resolutions to sweep.  Default: [0.2, 0.4, 0.6, 0.8, 1.0].
    best_resolution_override : float, optional
        Pin a specific resolution — bypasses all scoring criteria.
        Must be present in resolution_range.  Default: None.
    n_clusters_expected : int, optional
        Biological prior on the number of expected clusters (e.g. from
        previous work on the same tissue).  When provided, selects the
        resolution whose cluster count is closest to this target; uses
        the stability score as a tiebreaker.  Default: None.
    pca_key : str
        Key in adata.obsm to use for silhouette scoring.
    connectivities_key : str
        Key in adata.obsp for the Leiden adjacency matrix.
    random_state : int
        Reproducibility seed.
    inplace : bool
        If False (default) operate on a copy and leave the caller unchanged.

    Returns
    -------
    adata_out : AnnData
        AnnData with the following additions:
          - obs[f'leiden_{res}']      cluster labels for every resolution
          - obs['leiden']             labels at the selected resolution
          - uns['omicsage_cluster']   provenance record
    metrics : dict
        resolutions, n_clusters, n_clusters_delta, silhouette_scores,
        stability_scores, best_resolution, best_n_clusters,
        best_silhouette, selection_reason
    """
    _validate_inputs(adata, pca_key, connectivities_key, best_resolution_override,
                     resolution_range)

    if not inplace:
        adata = adata.copy()

    resolution_range = sorted(set(resolution_range))

    logger.info(
        "Running Leiden clustering — %d cells, %d resolutions: %s",
        adata.n_obs, len(resolution_range), resolution_range,
    )

    # ------------------------------------------------------------------
    # Sweep all resolutions
    # ------------------------------------------------------------------
    n_clusters:       dict[float, int]   = {}
    silhouette_scores: dict[float, float] = {}

    pca_matrix    = adata.obsm[pca_key]
    subsample_idx = _subsample_idx(adata.n_obs, _SILHOUETTE_MAX_CELLS, random_state)

    if subsample_idx is not None:
        logger.info(
            "Silhouette scoring will use a subsample of %d / %d cells",
            len(subsample_idx), adata.n_obs,
        )

    for res in resolution_range:
        obs_key = _res_key(res)
        sc.tl.leiden(
            adata,
            resolution=res,
            adjacency=adata.obsp[connectivities_key],
            key_added=obs_key,
            random_state=random_state,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        labels   = adata.obs[obs_key].values
        n_unique = len(np.unique(labels))
        n_clusters[res]        = n_unique
        silhouette_scores[res] = _silhouette(pca_matrix, labels, subsample_idx, n_unique)

        logger.info(
            "  resolution=%.2f → %d clusters, silhouette=%.4f",
            res, n_unique, silhouette_scores[res],
        )

    # ------------------------------------------------------------------
    # Derived metrics
    # ------------------------------------------------------------------
    # Delta: how many NEW clusters each resolution step adds vs the previous
    n_clusters_delta: dict[float, int] = {}
    for i, res in enumerate(resolution_range):
        if i == 0:
            n_clusters_delta[res] = 0          # no previous step
        else:
            prev = resolution_range[i - 1]
            n_clusters_delta[res] = n_clusters[res] - n_clusters[prev]

    # Stability score: 1 / (1 + delta)  →  1.0 at plateau, lower when jumping
    # Use the NEXT step's delta so that a resolution is rewarded for being
    # the last one before a plateau begins.
    stability_scores: dict[float, float] = {}
    for i, res in enumerate(resolution_range):
        if i < len(resolution_range) - 1:
            next_res   = resolution_range[i + 1]
            next_delta = abs(n_clusters_delta[next_res])
        else:
            next_delta = 0   # last resolution — no further jump to penalise
        stability_scores[res] = 1.0 / (1.0 + next_delta)

    logger.info("Cluster count deltas  : %s", n_clusters_delta)
    logger.info("Stability scores      : %s",
                {r: round(s, 3) for r, s in stability_scores.items()})

    # ------------------------------------------------------------------
    # Select best resolution
    # ------------------------------------------------------------------
    best_res, selection_reason = _select_resolution(
        resolution_range=resolution_range,
        n_clusters=n_clusters,
        n_clusters_delta=n_clusters_delta,
        silhouette_scores=silhouette_scores,
        stability_scores=stability_scores,
        best_resolution_override=best_resolution_override,
        n_clusters_expected=n_clusters_expected,
    )

    adata.obs["leiden"] = adata.obs[_res_key(best_res)].copy()

    logger.info(
        "Selected resolution: %.2f  (%d clusters, silhouette=%.4f, reason='%s')",
        best_res, n_clusters[best_res], silhouette_scores[best_res], selection_reason,
    )

    # ------------------------------------------------------------------
    # Build metrics dict (float keys for caller) + uns (string keys for h5ad)
    # ------------------------------------------------------------------
    metrics = {
        "resolutions":       resolution_range,
        "n_clusters":        n_clusters,
        "n_clusters_delta":  n_clusters_delta,
        "silhouette_scores": silhouette_scores,
        "stability_scores":  stability_scores,
        "best_resolution":   best_res,
        "best_n_clusters":   n_clusters[best_res],
        "best_silhouette":   silhouette_scores[best_res],
        "best_stability":    stability_scores[best_res],
        "selection_reason":  selection_reason,
        "n_clusters_expected": n_clusters_expected,
    }

    # HDF5 requires string dict keys
    def _str_keys(d: dict) -> dict:
        return {str(k): v for k, v in d.items()}

    adata.uns["omicsage_cluster"] = {
        "resolutions":         resolution_range,
        "n_clusters":          _str_keys(n_clusters),
        "n_clusters_delta":    _str_keys(n_clusters_delta),
        "silhouette_scores":   _str_keys(silhouette_scores),
        "stability_scores":    _str_keys(stability_scores),
        "best_resolution":     best_res,
        "best_n_clusters":     n_clusters[best_res],
        "best_silhouette":     silhouette_scores[best_res],
        "best_stability":      stability_scores[best_res],
        "selection_reason":    selection_reason,
        "n_clusters_expected": n_clusters_expected,
        "pca_key":             pca_key,
        "connectivities_key":  connectivities_key,
        "random_state":        random_state,
        "scanpy_version":      importlib.metadata.version("scanpy"),
        "omicsage_module":     "pipeline.modules.qc.cluster",
        "omicsage_version":    "0.1.0",
        "timestamp":           datetime.now(timezone.utc).isoformat(),
    }

    logger.info(
        "Clustering complete — resolution=%.2f, %d clusters, reason='%s'",
        best_res, n_clusters[best_res], selection_reason,
    )

    return adata, metrics


# ---------------------------------------------------------------------------
# Resolution selection logic
# ---------------------------------------------------------------------------

def _select_resolution(
    resolution_range:       list[float],
    n_clusters:             dict[float, int],
    n_clusters_delta:       dict[float, int],
    silhouette_scores:      dict[float, float],
    stability_scores:       dict[float, float],
    best_resolution_override: Optional[float],
    n_clusters_expected:    Optional[int],
) -> tuple[float, str]:
    """
    Select the best resolution using the priority order:
      1. best_resolution_override  (manual pin)
      2. n_clusters_expected       (closest count, stability tiebreaker)
      3. stability plateau         (silhouette tiebreaker)
      4. silhouette argmax         (fallback)

    Returns (best_res, selection_reason).
    """

    # Priority 1 — manual override
    if best_resolution_override is not None:
        return best_resolution_override, "override"

    # Priority 2 — biological prior (n_clusters_expected)
    if n_clusters_expected is not None:
        return _select_by_expected_count(
            resolution_range, n_clusters, stability_scores, n_clusters_expected
        )

    # Priority 3 — stability plateau
    return _select_by_stability(resolution_range, stability_scores, silhouette_scores)


def _select_by_expected_count(
    resolution_range:    list[float],
    n_clusters:          dict[float, int],
    stability_scores:    dict[float, float],
    n_clusters_expected: int,
) -> tuple[float, str]:
    """
    Pick the resolution whose cluster count is closest to n_clusters_expected.
    Break ties by highest stability score (i.e. prefer the plateau edge).
    """
    min_diff = min(abs(n_clusters[r] - n_clusters_expected) for r in resolution_range)
    candidates = [r for r in resolution_range
                  if abs(n_clusters[r] - n_clusters_expected) == min_diff]

    best_res = max(candidates, key=lambda r: stability_scores[r])

    logger.info(
        "Expected count selection: target=%d, candidates=%s (diff=%d), selected=%.2f",
        n_clusters_expected,
        {r: n_clusters[r] for r in candidates},
        min_diff,
        best_res,
    )
    return best_res, "expected_count"


def _select_by_stability(
    resolution_range:  list[float],
    stability_scores:  dict[float, float],
    silhouette_scores: dict[float, float],
) -> tuple[float, str]:
    """
    Pick the resolution at the elbow of the delta curve — the resolution
    where the cluster count makes its largest jump, after which further
    increases in resolution produce only small additions.

    This corresponds to the biological intuition: the resolution just after
    the biggest structural split captures the main cell populations without
    over-fragmenting into incremental sub-clusters.

    Example (from the session that motivated this fix):
      res=0.2 →  9 clusters  delta=0
      res=0.4 → 10 clusters  delta=+1
      res=0.6 → 17 clusters  delta=+7   ← largest jump → select 0.6
      res=0.8 → 20 clusters  delta=+3
      res=1.0 → 23 clusters  delta=+3
      res=1.2 → 26 clusters  delta=+3
      res=1.5 → 28 clusters  delta=+2   ← old (wrong) logic selected this

    Algorithm:
      1. Compute delta drop between consecutive resolutions:
         delta_drop[i] = delta[i] - delta[i+1]  (how much the jump shrinks)
      2. The resolution with the largest delta (most new clusters added) is
         the structural "unlock" point — select it.
      3. Use silhouette as tiebreaker if multiple resolutions share the max delta.
      4. Fall back to silhouette argmax if all deltas are equal (flat curve).
    """
    if len(resolution_range) < 2:
        best_res = resolution_range[0]
        return best_res, "silhouette"

    # deltas for resolutions after the first (first delta is always 0)
    deltas = {r: stability_scores[r] for r in resolution_range}  # kept for compat
    # Re-derive raw deltas from stability scores: stability = 1/(1+next_delta)
    # → next_delta = 1/stability - 1
    # Instead, recalculate directly from stability_scores which encode next_delta.
    # But stability_scores are already computed from n_clusters_delta in the caller.
    # We don't have n_clusters_delta here — use stability to reconstruct next_delta:
    #   stability[r] = 1 / (1 + |delta[next_r]|)
    #   |delta[next_r]| = 1/stability[r] - 1
    reconstructed_next_delta = {
        r: round(1.0 / stability_scores[r] - 1.0)
        for r in resolution_range
    }

    # The resolution whose NEXT step has the largest drop is the elbow:
    # we want the res where reconstructed_next_delta is largest (the big jump
    # happens going INTO the next resolution, so we pick the resolution BEFORE
    # the jump — that is, the resolution at which the next step adds the most).
    # But we want the resolution AFTER the big jump (the one that captured it).
    # The delta of res[i] = clusters added going from res[i-1] to res[i].
    # Largest delta means res[i] is where the unlock happened — select res[i].

    # next_delta[r] tells us: if I go from r to next_r, how many clusters do I add?
    # We want the r for which next_delta is largest — meaning going TO next_r was
    # the big jump — so select next_r (the resolution that made the big addition).
    max_next_delta = max(reconstructed_next_delta[r] for r in resolution_range)
    all_equal = max_next_delta == 0 or len(set(reconstructed_next_delta.values())) == 1

    if all_equal:
        best_res = max(silhouette_scores, key=lambda r: silhouette_scores[r])
        logger.info(
            "Delta curve is flat — falling back to silhouette argmax: %.2f", best_res
        )
        return best_res, "silhouette"

    # Find the resolution(s) with the largest next_delta
    candidates_before = [r for r in resolution_range
                         if reconstructed_next_delta[r] == max_next_delta]

    # Select the NEXT resolution after the biggest jump — that is where the
    # structural unlock was captured
    best_candidates = []
    for r in candidates_before:
        idx = resolution_range.index(r)
        if idx + 1 < len(resolution_range):
            best_candidates.append(resolution_range[idx + 1])
        else:
            best_candidates.append(r)  # edge case: last resolution

    # Tiebreak by silhouette
    best_res = max(best_candidates, key=lambda r: silhouette_scores[r])

    logger.info(
        "Delta-elbow selection: largest jump before=%s (delta=%d), "
        "selected=%s (silhouette=%.4f)",
        candidates_before, max_next_delta, best_res, silhouette_scores[best_res],
    )
    return best_res, "stability_plateau"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _res_key(res: float) -> str:
    """Convert a float resolution to a clean obs column name."""
    return f"leiden_{res:.4g}"


def _validate_inputs(
    adata: AnnData,
    pca_key: str,
    connectivities_key: str,
    best_resolution_override: Optional[float],
    resolution_range: list[float],
) -> None:
    if pca_key not in adata.obsm:
        raise KeyError(
            f"adata.obsm['{pca_key}'] not found. Run reduce() before cluster()."
        )
    if connectivities_key not in adata.obsp:
        raise KeyError(
            f"adata.obsp['{connectivities_key}'] not found. Run reduce() before cluster()."
        )
    if adata.n_obs < 2:
        raise ValueError("AnnData must contain at least 2 cells.")
    if best_resolution_override is not None:
        res_set = sorted(set(resolution_range))
        if best_resolution_override not in res_set:
            raise ValueError(
                f"best_resolution_override={best_resolution_override} is not in "
                f"resolution_range={res_set}. Add it to resolution_range or choose "
                "a value from that list."
            )


def _subsample_idx(
    n_obs: int,
    max_cells: int,
    random_state: int,
) -> Optional[np.ndarray]:
    if n_obs <= max_cells:
        return None
    rng = np.random.default_rng(random_state)
    return rng.choice(n_obs, size=max_cells, replace=False)


def _silhouette(
    pca_matrix: np.ndarray,
    labels: np.ndarray,
    subsample_idx: Optional[np.ndarray],
    n_unique: int,
) -> float:
    if n_unique < 2:
        return -1.0
    if subsample_idx is not None:
        X = pca_matrix[subsample_idx]
        y = labels[subsample_idx]
        if len(np.unique(y)) < 2:
            X, y = pca_matrix, labels
    else:
        X, y = pca_matrix, labels
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return float(silhouette_score(X, y, metric="euclidean"))
    except Exception:
        return -1.0
