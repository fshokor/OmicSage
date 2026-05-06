"""
pipeline/modules/qc/cluster.py
OmicSage — Leiden clustering module

Input  : reduced AnnData (output of reduce.py)
         Required: obsm['X_pca'], obsp['connectivities']
Output : (adata_clustered, metrics_dict)

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
_SILHOUETTE_MAX_CELLS = 10_000   # subsample threshold for silhouette computation


def cluster(
    adata: AnnData,
    resolution_range: list[float] = DEFAULT_RESOLUTIONS,
    best_resolution_override: Optional[float] = None,
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
    pca_key : str
        Key in adata.obsm to use as feature matrix for silhouette scoring.
    best_resolution_override : float, optional
        If provided, skip silhouette-based auto-selection and use this
        resolution as the best.  The resolution must be present in
        resolution_range — add it automatically if missing.
        Use this when you want more clusters than silhouette would select
        (e.g. to over-cluster before annotation).  Default: None (auto).
    connectivities_key : str
        Key in adata.obsp to use as the pre-built adjacency matrix.
    random_state : int
        Random seed for reproducibility.
    inplace : bool
        If False (default) operate on a copy and leave the caller unchanged.

    Returns
    -------
    adata_out : AnnData
        AnnData with the following additions:
          - obs[f'leiden_{res}']          cluster labels for every resolution
          - obs['leiden']                 labels at the auto-selected resolution
          - uns['omicsage_cluster']       provenance record
    metrics : dict
        {
          'resolutions'        : list of floats tested,
          'n_clusters'         : {res: int},
          'silhouette_scores'  : {res: float},
          'best_resolution'    : float,
          'best_silhouette'    : float,
          'best_n_clusters'    : int,
        }
    """
    _validate_inputs(adata, pca_key, connectivities_key)

    if not inplace:
        adata = adata.copy()

    resolution_range = sorted(set(resolution_range))   # deduplicate + sort

    logger.info(
        "Running Leiden clustering — %d cells, %d resolutions: %s",
        adata.n_obs,
        len(resolution_range),
        resolution_range,
    )

    n_clusters: dict[float, int] = {}
    silhouette_scores: dict[float, float] = {}

    pca_matrix = adata.obsm[pca_key]
    subsample_idx = _subsample_idx(adata.n_obs, _SILHOUETTE_MAX_CELLS, random_state)

    if subsample_idx is not None:
        logger.info(
            "Silhouette scoring will use a subsample of %d / %d cells",
            len(subsample_idx),
            adata.n_obs,
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

        labels = adata.obs[obs_key].values
        n_unique = len(np.unique(labels))
        n_clusters[res] = n_unique

        sil = _silhouette(pca_matrix, labels, subsample_idx, n_unique)
        silhouette_scores[res] = sil

        logger.info(
            "  resolution=%.2f → %d clusters, silhouette=%.4f",
            res,
            n_unique,
            sil,
        )

    # ------------------------------------------------------------------
    # Select best resolution — override or silhouette auto-select
    # ------------------------------------------------------------------
    if best_resolution_override is not None:
        if best_resolution_override not in silhouette_scores:
            raise ValueError(
                f"best_resolution_override={best_resolution_override} is not in "
                f"resolution_range={resolution_range}. "
                "Add it to resolution_range or choose a value from that list."
            )
        best_res = best_resolution_override
        logger.info(
            "Resolution overridden by user: %.2f (%d clusters, silhouette=%.4f)",
            best_res,
            n_clusters[best_res],
            silhouette_scores[best_res],
        )
    else:
        best_res = max(silhouette_scores, key=lambda r: silhouette_scores[r])
        logger.info(
            "Best resolution auto-selected: %.2f (%d clusters, silhouette=%.4f)",
            best_res,
            n_clusters[best_res],
            silhouette_scores[best_res],
        )

    adata.obs["leiden"] = adata.obs[_res_key(best_res)].copy()

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------
    metrics = {
        "resolutions": resolution_range,
        "n_clusters": n_clusters,
        "silhouette_scores": silhouette_scores,
        "best_resolution": best_res,
        "best_silhouette": silhouette_scores[best_res],
        "best_n_clusters": n_clusters[best_res],
    }

    # HDF5 requires string keys — convert float resolution keys before storing in uns
    adata.uns["omicsage_cluster"] = {
        "resolutions": resolution_range,
        "n_clusters": {str(r): v for r, v in n_clusters.items()},
        "silhouette_scores": {str(r): v for r, v in silhouette_scores.items()},
        "best_resolution": best_res,
        "best_silhouette": silhouette_scores[best_res],
        "best_n_clusters": n_clusters[best_res],
        "resolution_selection": "override" if best_resolution_override is not None else "silhouette",
        "pca_key": pca_key,
        "connectivities_key": connectivities_key,
        "random_state": random_state,
        "scanpy_version": importlib.metadata.version("scanpy"),
        "omicsage_module": "pipeline.modules.qc.cluster",
        "omicsage_version": "0.1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    logger.info(
        "Clustering complete — best resolution=%.2f, %d clusters",
        best_res,
        n_clusters[best_res],
    )

    return adata, metrics


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _res_key(res: float) -> str:
    """Convert a float resolution to a clean obs column name."""
    # Use up to 4 significant digits, strip trailing zeros
    return f"leiden_{res:.4g}"


def _validate_inputs(
    adata: AnnData,
    pca_key: str,
    connectivities_key: str,
) -> None:
    if pca_key not in adata.obsm:
        raise KeyError(
            f"adata.obsm['{pca_key}'] not found. "
            "Run reduce() before cluster()."
        )
    if connectivities_key not in adata.obsp:
        raise KeyError(
            f"adata.obsp['{connectivities_key}'] not found. "
            "Run reduce() before cluster()."
        )
    if adata.n_obs < 2:
        raise ValueError("AnnData must contain at least 2 cells.")


def _subsample_idx(
    n_obs: int,
    max_cells: int,
    random_state: int,
) -> Optional[np.ndarray]:
    """Return a random subsample index array, or None if no subsampling needed."""
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
    """
    Compute silhouette score.  Returns -1.0 if fewer than 2 clusters
    (silhouette is undefined) or if computation raises any warning/error.
    """
    if n_unique < 2:
        return -1.0

    if subsample_idx is not None:
        X = pca_matrix[subsample_idx]
        y = labels[subsample_idx]
        # Ensure subsampled labels still have ≥2 unique values
        if len(np.unique(y)) < 2:
            X = pca_matrix
            y = labels
    else:
        X = pca_matrix
        y = labels

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return float(silhouette_score(X, y, metric="euclidean"))
    except Exception:
        return -1.0

