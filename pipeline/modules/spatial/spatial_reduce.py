"""
spatial_reduce.py — OmicSage Phase 7, Session 2
Normalization, HVG selection, PCA, and spatial neighbours for Visium spots.

Input:  AnnData from spatial_qc (has obsm["spatial"], qc_pass column)
Output: AnnData with normalized X, obsm["X_pca"], spatial graph in
        obsp["spatial_connectivities"], provenance in
        uns["omicsage_spatial_reduce"]

Raw counts are preserved in layers["counts"] throughout.

Notes
-----
- The squidpy benchmark dataset (sq.datasets.visium_hne_adata()) is
  pre-processed (already normalized + log1p + clustered).  We detect
  this by checking layers["counts"]: if absent AND the ingest provenance
  says spatial_type="benchmark", normalization is skipped.
  The safest general rule: if layers["counts"] is present, normalize
  from there.  If absent we normalize adata.X in-place (raw Visium path).
- coord_type=None lets squidpy auto-detect "grid" when uns["spatial"]
  is present, which is correct for standard Visium.  Pass coord_type
  explicitly only to override for non-grid data.
- flavor="seurat" for HVG is safe for RNA (no TF-IDF overflow risk).
"""

from __future__ import annotations

from datetime import datetime
from typing import Optional

import anndata as ad
import numpy as np
import scanpy as sc

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spatial_reduce(
    adata: ad.AnnData,
    n_top_genes: int = 3000,
    n_comps: int = 50,
    n_neighbors: int = 6,
    coord_type: Optional[str] = None,
    normalize_total: bool = True,
    target_sum: float = 1e4,
    log1p: bool = True,
    flavor: str = "seurat",
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Normalize, select HVGs, compute PCA and spatial neighbours graph.

    Parameters
    ----------
    adata
        AnnData produced by :func:`spatial_qc`.  Must contain
        ``obsm["spatial"]``.
    n_top_genes
        Number of highly variable genes to select.
        Default 3000 is appropriate for Visium (~33k genes).
    n_comps
        Number of PCA components to compute.
    n_neighbors
        Number of spatial neighbours for the neighbours graph.
        Default 6 matches the Visium hexagonal grid topology.
        Only used when ``coord_type`` is ``"grid"`` or ``None``
        (auto-detected as grid when ``uns["spatial"]`` is present).
    coord_type
        Coordinate type passed to ``squidpy.gr.spatial_neighbors``.
        ``None`` (default) — auto-detect; uses ``"grid"`` when
        ``uns["spatial"]`` is present (standard Visium).
        ``"grid"`` — force hexagonal grid connectivity.
        ``"generic"`` — generic (non-grid) connectivity.
    normalize_total
        Whether to run ``sc.pp.normalize_total``.
        Set to ``False`` if data is already normalized.
    target_sum
        Normalization target (CPM-like). Default ``1e4``.
    log1p
        Whether to apply ``sc.pp.log1p`` after normalization.
    flavor
        HVG selection flavor: ``"seurat"`` (default) or
        ``"cell_ranger"``.  ``"seurat"`` is safe for RNA data.
    inplace
        If ``False`` (default), operate on a copy and return it.
        If ``True``, modify *adata* in place.

    Returns
    -------
    adata
        AnnData with:

        - ``X`` — normalized + log1p counts
        - ``layers["counts"]`` — raw integer counts (preserved)
        - ``var["highly_variable"]`` — HVG boolean mask
        - ``obsm["X_pca"]`` — PCA embedding
        - ``uns["pca"]["variance_ratio"]`` — explained variance per PC
        - ``obsp["spatial_connectivities"]`` — spatial adjacency
        - ``obsp["spatial_distances"]`` — spatial distances
        - ``uns["spatial_neighbors"]`` — spatial graph metadata
        - ``uns["omicsage_spatial_reduce"]`` — provenance

    params
        Provenance dictionary (same as ``uns["omicsage_spatial_reduce"]``).

    Raises
    ------
    ImportError
        If squidpy is not installed.
    TypeError
        If *adata* is not an AnnData object.
    ValueError
        If ``obsm["spatial"]`` is missing.
    """
    _check_squidpy()
    _validate_input(adata)

    adata = adata if inplace else adata.copy()

    # ------------------------------------------------------------------ #
    # 1.  Preserve raw counts in layers["counts"]
    # ------------------------------------------------------------------ #
    already_normalized = False

    if "counts" in adata.layers:
        # Raw counts already saved (visium path via spatial_ingest, or
        # caller preserved them).  Restore raw counts to X before normalizing.
        adata.X = adata.layers["counts"].copy()
    else:
        # Check whether this is a pre-processed benchmark dataset.
        ingest_meta = adata.uns.get("omicsage_spatial_ingest", {})
        spatial_type = ingest_meta.get("spatial_type", "")
        if spatial_type == "benchmark":
            # Benchmark data is already normalized — do not re-normalize.
            already_normalized = True
        # Save whatever is in X as "counts" for downstream reference.
        # For benchmark data this is normalized values, not raw integers.
        adata.layers["counts"] = adata.X.copy()

    # ------------------------------------------------------------------ #
    # 2.  Normalize + log1p
    # ------------------------------------------------------------------ #
    if not already_normalized and normalize_total:
        sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

    if not already_normalized and log1p:
        sc.pp.log1p(adata)

    # ------------------------------------------------------------------ #
    # 3.  Highly variable gene selection
    # ------------------------------------------------------------------ #
    n_top_genes_actual = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes_actual,
        flavor=flavor,
        inplace=True,
    )
    n_hvgs = int(adata.var["highly_variable"].sum())

    # ------------------------------------------------------------------ #
    # 4.  PCA (on HVG subset)
    # ------------------------------------------------------------------ #
    n_comps_actual = min(n_comps, n_hvgs - 1, adata.n_obs - 1)
    sc.pp.pca(
        adata,
        n_comps=n_comps_actual,
        use_highly_variable=True,
        svd_solver="arpack",
    )
    variance_ratio = adata.uns["pca"]["variance_ratio"].tolist()

    # ------------------------------------------------------------------ #
    # 5.  Spatial neighbours graph
    # ------------------------------------------------------------------ #
    sq.gr.spatial_neighbors(
        adata,
        n_neighs=n_neighbors,
        coord_type=coord_type,
        key_added="spatial",
    )

    # Connectivity stats
    conn = adata.obsp["spatial_connectivities"]
    n_edges = int(conn.nnz // 2)
    mean_neighbors = float(np.array(conn.sum(axis=1)).mean())

    # ------------------------------------------------------------------ #
    # 6.  Provenance
    # ------------------------------------------------------------------ #
    params = {
        "module": "spatial_reduce",
        "timestamp": datetime.now().isoformat(),
        "params": {
            "n_top_genes": n_top_genes,
            "n_comps": n_comps,
            "n_neighbors": n_neighbors,
            "coord_type": coord_type,
            "normalize_total": normalize_total,
            "target_sum": target_sum,
            "log1p": log1p,
            "flavor": flavor,
            "skipped_normalization": already_normalized,
        },
        "outputs": {
            "n_hvgs": n_hvgs,
            "n_comps_computed": n_comps_actual,
            "pca_variance_ratio_top10": variance_ratio[:10],
            "pca_cumulative_variance_top10": float(
                np.cumsum(variance_ratio[:10])[-1]
            ),
            "spatial_graph_n_edges": n_edges,
            "spatial_graph_mean_neighbors": round(mean_neighbors, 2),
        },
    }
    adata.uns["omicsage_spatial_reduce"] = params
    return adata, params


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _check_squidpy() -> None:
    if not _SQUIDPY_AVAILABLE:
        raise ImportError(
            "squidpy is required for spatial analysis. "
            "Install with: pip install squidpy"
        )


def _validate_input(adata: ad.AnnData) -> None:
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"Expected AnnData, got {type(adata).__name__}")
    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] is missing. "
            "Run spatial_ingest() and spatial_qc() first."
        )
    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 observations (spots). Cannot run reduce.")
    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 variables (genes). Cannot run reduce.")
