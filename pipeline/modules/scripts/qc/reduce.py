"""
OmicSage — Dimensionality Reduction Module
pipeline/modules/qc/reduce.py

Runs PCA, neighbor graph construction, UMAP (and optional t-SNE) on a
normalized AnnData (output of normalize.py).

The number of PCs used for the neighbor graph is chosen automatically by
default using elbow detection on the variance-explained curve (kneed).
The user can override this with n_pcs=<int>.

Input : normalized AnnData (output of normalize.py)
        - log1p values in .X
        - HVGs flagged in .var['highly_variable']
Output: reduced AnnData + metrics dict

Steps
-----
1. PCA on HVG subset (n_comps=50 by default)
2. Auto-select number of PCs to use for neighbors (elbow / variance / fixed)
3. Neighbor graph (sc.pp.neighbors)
4. UMAP (sc.tl.umap) — always on
5. t-SNE (sc.tl.tsne) — optional, default off
6. Store all params + provenance in adata.uns['omicsage_reduce']

Usage
-----
    from pipeline.modules.qc.reduce import reduce

    adata_reduced, metrics = reduce(adata_norm)

    # Manual PC override:
    adata_reduced, metrics = reduce(adata_norm, n_pcs=20)

    # With t-SNE:
    adata_reduced, metrics = reduce(adata_norm, run_tsne=True)
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import Literal, Optional

import numpy as np
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def reduce(
    adata: AnnData,
    n_comps: int = 50,
    n_pcs: Optional[int] = None,
    n_pcs_method: Literal["elbow", "variance", "fixed"] = "elbow",
    variance_threshold: float = 0.85,
    n_neighbors: int = 15,
    run_tsne: bool = False,
    inplace: bool = False,
    random_state: int = 0,
) -> tuple[AnnData, dict]:
    """
    Run PCA, neighbor graph, UMAP (and optionally t-SNE) on a normalized AnnData.

    Parameters
    ----------
    adata : AnnData
        Normalized AnnData with log1p values in .X and HVGs flagged in
        .var['highly_variable'].  Typically the output of normalize().
    n_comps : int
        Number of PCA components to compute.  Default 50.
    n_pcs : int, optional
        Number of PCs to pass to the neighbor graph.  If None (default),
        the value is chosen automatically via n_pcs_method.
        Pass an explicit integer to override auto-selection entirely.
    n_pcs_method : {'elbow', 'variance', 'fixed'}
        Strategy used when n_pcs is None:
        - 'elbow'    : kneed elbow detection on the variance-explained curve.
                       Falls back to 'variance' if kneed is unavailable or
                       fails to find a clear elbow.
        - 'variance' : keep PCs until cumulative variance explained >=
                       variance_threshold.
        - 'fixed'    : always use min(30, n_comps).
        Default 'elbow'.
    variance_threshold : float
        Cumulative variance explained target for n_pcs_method='variance'.
        Default 0.85 (85 %).
    n_neighbors : int
        Number of neighbors for the kNN graph.  Default 15.
    run_tsne : bool
        If True, also compute t-SNE and store in obsm['X_tsne'].
        Default False (t-SNE is slow on large datasets).
    inplace : bool
        If False (default) work on a copy so the caller's object is not mutated.
    random_state : int
        Reproducibility seed passed to all stochastic steps.

    Returns
    -------
    adata_out : AnnData
        AnnData with embeddings populated:
        - obsm['X_pca']    — PCA embedding (n_cells × n_comps)
        - obsm['X_umap']   — UMAP embedding (n_cells × 2)
        - obsm['X_tsne']   — t-SNE embedding (n_cells × 2) if run_tsne=True
        - obsp['connectivities'], obsp['distances'] — kNN graph
        - uns['omicsage_reduce'] — full provenance record
    metrics : dict
        Summary statistics for logging / downstream QA.
    """
    _validate_input(adata)

    # Work on a copy unless caller opts in to inplace mutation
    adata_out = adata if inplace else adata.copy()

    # ------------------------------------------------------------------
    # Step 1 — PCA on HVG subset
    # ------------------------------------------------------------------
    use_hvg = "highly_variable" in adata_out.var.columns and adata_out.var["highly_variable"].sum() > 0
    if not use_hvg:
        logger.warning(
            "No HVGs found in adata.var['highly_variable']. "
            "Running PCA on all genes. Consider running normalize() first."
        )

    logger.info(
        "Running PCA — n_comps=%d, use_highly_variable=%s",
        n_comps,
        use_hvg,
    )
    sc.tl.pca(
        adata_out,
        n_comps=n_comps,
        use_highly_variable=use_hvg,
        svd_solver="arpack",
        random_state=random_state,
    )

    # ------------------------------------------------------------------
    # Step 2 — Auto-select number of PCs for neighbor graph
    # ------------------------------------------------------------------
    variance_ratio = adata_out.uns["pca"]["variance_ratio"]   # shape (n_comps,)

    if n_pcs is not None:
        # User override — skip auto-selection entirely
        n_pcs_used = int(np.clip(n_pcs, 1, n_comps))
        pc_selection_method = "manual"
        logger.info("Using manually specified n_pcs=%d", n_pcs_used)
    else:
        n_pcs_used, pc_selection_method = _select_n_pcs(
            variance_ratio=variance_ratio,
            method=n_pcs_method,
            variance_threshold=variance_threshold,
            n_comps=n_comps,
        )

    logger.info(
        "PC selection — method='%s', n_pcs_used=%d "
        "(cumulative variance explained: %.1f%%)",
        pc_selection_method,
        n_pcs_used,
        float(variance_ratio[:n_pcs_used].sum()) * 100,
    )

    # ------------------------------------------------------------------
    # Step 3 — Neighbor graph
    # ------------------------------------------------------------------
    logger.info(
        "Computing neighbor graph — n_neighbors=%d, n_pcs=%d",
        n_neighbors,
        n_pcs_used,
    )
    sc.pp.neighbors(
        adata_out,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs_used,
        random_state=random_state,
    )

    # ------------------------------------------------------------------
    # Step 4 — UMAP
    # ------------------------------------------------------------------
    logger.info("Running UMAP")
    sc.tl.umap(adata_out, random_state=random_state)

    # ------------------------------------------------------------------
    # Step 5 — t-SNE (optional)
    # ------------------------------------------------------------------
    if run_tsne:
        logger.info("Running t-SNE (n_pcs=%d)", n_pcs_used)
        sc.tl.tsne(
            adata_out,
            n_pcs=n_pcs_used,
            random_state=random_state,
            use_rep="X_pca",
        )
    else:
        logger.debug("t-SNE skipped (run_tsne=False)")

    # ------------------------------------------------------------------
    # Step 6 — Build metrics dict
    # ------------------------------------------------------------------
    metrics = {
        "n_cells": int(adata_out.n_obs),
        "n_genes": int(adata_out.n_vars),
        "n_hvg_used": int(adata_out.var["highly_variable"].sum()) if use_hvg else int(adata_out.n_vars),
        "n_comps_computed": n_comps,
        "n_pcs_used": n_pcs_used,
        "pc_selection_method": pc_selection_method,
        "n_neighbors": n_neighbors,
        "variance_explained_per_pc": variance_ratio.tolist(),
        "cumulative_variance_explained": float(variance_ratio[:n_pcs_used].sum()),
        "run_tsne": run_tsne,
        "embeddings_computed": ["X_pca", "X_umap"] + (["X_tsne"] if run_tsne else []),
    }

    # ------------------------------------------------------------------
    # Step 7 — Store provenance in uns
    # ------------------------------------------------------------------
    adata_out.uns["omicsage_reduce"] = {
        "n_comps": n_comps,
        "n_pcs_used": n_pcs_used,
        "pc_selection_method": pc_selection_method,
        "variance_threshold": variance_threshold,
        "n_neighbors": n_neighbors,
        "run_tsne": run_tsne,
        "use_highly_variable": use_hvg,
        "n_hvg_used": metrics["n_hvg_used"],
        "variance_explained_per_pc": variance_ratio.tolist(),
        "cumulative_variance_explained_by_selected_pcs": metrics["cumulative_variance_explained"],
        "random_state": random_state,
        "scanpy_version": sc.__version__,
        "omicsage_module": "pipeline.modules.qc.reduce",
        "omicsage_version": "0.1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    logger.info(
        "Dimensionality reduction complete — "
        "%d cells, %d PCs used, UMAP computed%s",
        adata_out.n_obs,
        n_pcs_used,
        ", t-SNE computed" if run_tsne else "",
    )

    return adata_out, metrics


# ---------------------------------------------------------------------------
# PC selection
# ---------------------------------------------------------------------------

def _select_n_pcs(
    variance_ratio: np.ndarray,
    method: str,
    variance_threshold: float,
    n_comps: int,
) -> tuple[int, str]:
    """
    Choose the number of PCs to use for the neighbor graph.

    Returns
    -------
    n_pcs : int
        Chosen number of PCs.
    method_used : str
        Label describing which method was actually used (may differ from
        *method* if a fallback was triggered).
    """
    if method == "elbow":
        n_pcs, method_used = _select_elbow(variance_ratio, n_comps)
        return n_pcs, method_used

    if method == "variance":
        return _select_variance(variance_ratio, variance_threshold, n_comps)

    # method == "fixed"
    return _select_fixed(n_comps)


def _select_elbow(
    variance_ratio: np.ndarray,
    n_comps: int,
) -> tuple[int, str]:
    """
    Use the kneed KneeLocator to find the elbow of the scree plot.

    Falls back to variance threshold (85%) if kneed is not installed or
    fails to find a clear elbow.
    """
    try:
        from kneed import KneeLocator

        x = np.arange(1, len(variance_ratio) + 1)
        y = variance_ratio

        kneedle = KneeLocator(
            x,
            y,
            S=1.0,               # sensitivity — higher = less sensitive
            curve="convex",
            direction="decreasing",
            interp_method="interp1d",
        )

        if kneedle.elbow is not None:
            n_pcs = int(kneedle.elbow)
            # Enforce reasonable bounds: at least 5 PCs, at most n_comps
            n_pcs = int(np.clip(n_pcs, 5, n_comps))
            logger.info("Elbow detected at PC %d", n_pcs)
            return n_pcs, "elbow"
        else:
            logger.warning(
                "kneed could not detect a clear elbow — "
                "falling back to variance threshold method."
            )
            return _select_variance(variance_ratio, 0.85, n_comps)

    except ImportError:
        logger.warning(
            "kneed package not found — falling back to variance threshold method. "
            "Install with: pip install kneed"
        )
        return _select_variance(variance_ratio, 0.85, n_comps)

    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "Elbow detection failed (%s) — falling back to variance threshold method.", exc
        )
        return _select_variance(variance_ratio, 0.85, n_comps)


def _select_variance(
    variance_ratio: np.ndarray,
    variance_threshold: float,
    n_comps: int,
) -> tuple[int, str]:
    """
    Keep PCs until cumulative variance explained >= variance_threshold.

    Always returns at least 5 PCs and at most n_comps PCs.
    """
    cumvar = np.cumsum(variance_ratio)
    above = np.where(cumvar >= variance_threshold)[0]

    if len(above) == 0:
        # Variance threshold not reached — use all computed PCs
        logger.warning(
            "Cumulative variance never reached %.0f%% across %d PCs — "
            "using all %d PCs.",
            variance_threshold * 100,
            n_comps,
            n_comps,
        )
        n_pcs = n_comps
    else:
        n_pcs = int(above[0]) + 1   # +1 because indices are 0-based

    n_pcs = int(np.clip(n_pcs, 5, n_comps))
    logger.info(
        "Variance threshold method: %d PCs explain %.1f%% of variance",
        n_pcs,
        float(cumvar[n_pcs - 1]) * 100,
    )
    return n_pcs, "variance"


def _select_fixed(n_comps: int) -> tuple[int, str]:
    """Return a fixed default of min(30, n_comps) PCs."""
    n_pcs = min(30, n_comps)
    logger.info("Fixed PC selection: n_pcs=%d", n_pcs)
    return n_pcs, "fixed"


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def _validate_input(adata: AnnData) -> None:
    """Raise informative errors for common input mistakes."""
    if not isinstance(adata, AnnData):
        raise TypeError(
            f"reduce() expects an AnnData object, got {type(adata).__name__}. "
            "Pass mdata['rna'] not the full MuData."
        )

    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 cells — nothing to reduce.")

    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 genes — nothing to reduce.")

    # Check .X looks like log-normalized data (not raw counts)
    # Raw counts would have many exact integers and a very different distribution
    x_sample = adata.X[:min(200, adata.n_obs), :min(200, adata.n_vars)]
    if hasattr(x_sample, "toarray"):
        x_sample = x_sample.toarray()
    x_sample = np.asarray(x_sample, dtype=float)

    all_integer_frac = np.mean(np.abs(x_sample - np.round(x_sample)) < 1e-6)
    if all_integer_frac > 0.95 and x_sample.max() > 50:
        logger.warning(
            "adata.X appears to contain raw counts (%.0f%% integer values, max=%.0f). "
            "reduce() expects log-normalized data. "
            "Did you run normalize() first?",
            all_integer_frac * 100,
            x_sample.max(),
        )

    logger.debug("Input validation passed — %d cells × %d genes", adata.n_obs, adata.n_vars)
