"""
harmony_correct.py — Batch correction via Harmony for OmicSage Phase 1.

Workflow
--------
1. Validate that PCA has been run (obsm['X_pca'] must exist).
2. Run Harmony on X_pca using the specified batch_key column in obs.
3. Store corrected embedding in obsm['X_pca_harmony'].
4. Recompute neighbor graph and UMAP on the corrected embedding.
5. Write provenance to uns['omicsage_harmony'].

Public API
----------
    harmony_correct(adata, batch_key="batch", ...) -> AnnData
"""

from __future__ import annotations

import logging
import time
from typing import Optional

import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def harmony_correct(
    adata,
    batch_key: str = "batch",
    n_pcs: int = 50,
    n_neighbors: int = 15,
    umap_min_dist: float = 0.3,
    random_state: int = 0,
    max_iter_harmony: int = 50,
    theta: float = 2.0,
    copy: bool = False,
) -> "AnnData":
    """Run Harmony batch correction and recompute graph + UMAP.

    Parameters
    ----------
    adata : AnnData
        Must have ``obsm['X_pca']`` and ``obs[batch_key]``.
    batch_key : str
        Column in ``adata.obs`` that encodes the batch variable.
    n_pcs : int
        Number of PCs to pass to Harmony (uses the first *n_pcs* dims of
        ``X_pca``; capped at the actual number of PCs available).
    n_neighbors : int
        *k* for the neighbor graph built on the corrected embedding.
    umap_min_dist : float
        ``min_dist`` parameter for UMAP.
    random_state : int
        Random seed for reproducibility.
    max_iter_harmony : int
        Maximum number of Harmony iterations.
    theta : float
        Harmony diversity penalty (higher = stronger correction).
    copy : bool
        If True, return a copy; otherwise modify in place.

    Returns
    -------
    AnnData
        With ``obsm['X_pca_harmony']``, updated neighbor graph, and
        ``obsm['X_umap']`` (recomputed on the corrected embedding).
    """
    adata = adata.copy() if copy else adata

    _validate_inputs(adata, batch_key)

    n_pcs_actual = min(n_pcs, adata.obsm["X_pca"].shape[1])
    if n_pcs_actual < n_pcs:
        logger.warning(
            "Requested %d PCs but X_pca only has %d — using %d.",
            n_pcs, adata.obsm["X_pca"].shape[1], n_pcs_actual,
        )

    t0 = time.time()

    # ------------------------------------------------------------------ #
    # 1. Harmony                                                           #
    # ------------------------------------------------------------------ #
    logger.info("Running Harmony on batch_key='%s' (%d PCs)…", batch_key, n_pcs_actual)
    _run_harmony(adata, batch_key=batch_key, n_pcs=n_pcs_actual,
                 max_iter=max_iter_harmony, theta=theta,
                 random_state=random_state)

    # ------------------------------------------------------------------ #
    # 2. Neighbor graph on corrected embedding                            #
    # ------------------------------------------------------------------ #
    logger.info("Computing neighbor graph on X_pca_harmony (k=%d)…", n_neighbors)
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep="X_pca_harmony",
        random_state=random_state,
        key_added="neighbors_harmony",
    )

    # ------------------------------------------------------------------ #
    # 3. Preserve existing UMAP, then compute new one on corrected graph  #
    # ------------------------------------------------------------------ #
    if "X_umap" in adata.obsm:
        adata.obsm["X_umap_precorrection"] = adata.obsm["X_umap"].copy()
        logger.info("Saved existing X_umap → X_umap_precorrection.")

    logger.info("Computing UMAP on corrected graph…")
    sc.tl.umap(
        adata,
        min_dist=umap_min_dist,
        random_state=random_state,
        neighbors_key="neighbors_harmony",
    )
    # Rename so callers can distinguish pre- vs post-correction UMAPs
    adata.obsm["X_umap_harmony"] = adata.obsm.pop("X_umap")

    elapsed = time.time() - t0

    # ------------------------------------------------------------------ #
    # 4. Provenance                                                        #
    # ------------------------------------------------------------------ #
    n_batches = int(adata.obs[batch_key].nunique())
    batch_sizes = (
        adata.obs[batch_key]
        .value_counts()
        .to_dict()
    )
    # JSON-safe: all keys to str
    batch_sizes = {str(k): int(v) for k, v in batch_sizes.items()}

    adata.uns["omicsage_harmony"] = {
        "batch_key": batch_key,
        "n_batches": n_batches,
        "batch_sizes": batch_sizes,
        "n_pcs": n_pcs_actual,
        "n_neighbors": n_neighbors,
        "umap_min_dist": umap_min_dist,
        "max_iter_harmony": max_iter_harmony,
        "theta": theta,
        "random_state": random_state,
        "elapsed_seconds": round(elapsed, 2),
        "embedding_key": "X_pca_harmony",
        "neighbors_key": "neighbors_harmony",
        "umap_key": "X_umap_harmony",
        "umap_precorrection_key": "X_umap_precorrection",
    }

    logger.info(
        "Harmony complete: %d batches, %.1fs elapsed.", n_batches, elapsed
    )
    return adata


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_inputs(adata, batch_key: str) -> None:
    """Raise informative errors for missing prerequisites."""
    if "X_pca" not in adata.obsm:
        raise ValueError(
            "obsm['X_pca'] not found. Run dimensionality reduction "
            "(reduce.py) before batch correction."
        )
    if batch_key not in adata.obs.columns:
        raise KeyError(
            f"batch_key='{batch_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        raise ValueError(
            f"batch_key='{batch_key}' has only {n_batches} unique value(s). "
            "Harmony requires at least 2 batches."
        )


def _run_harmony(adata, batch_key: str, n_pcs: int,
                 max_iter: int, theta: float, random_state: int) -> None:
    """Run harmonypy and store result in obsm['X_pca_harmony']."""
    try:
        import harmonypy as hm
    except ImportError as exc:
        raise ImportError(
            "harmonypy is not installed. Run:\n"
            "    pip install harmonypy\n"
            "then retry."
        ) from exc

    pca_matrix = adata.obsm["X_pca"][:, :n_pcs].copy()
    meta = adata.obs[[batch_key]].copy()

    ho = hm.run_harmony(
        pca_matrix,
        meta,
        batch_key,
        max_iter_harmony=max_iter,
        theta=theta,
        random_state=random_state,
    )

    # ho.Z_corr is (n_cells, n_pcs) in this version of harmonypy —
    # store directly; no transpose needed.
    adata.obsm["X_pca_harmony"] = np.array(ho.Z_corr).astype(np.float32)
