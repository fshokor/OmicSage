"""
spatial_cluster.py — OmicSage Phase 7, Session 3
Leiden clustering and spatially variable gene detection for Visium spots.

Input:  AnnData from spatial_reduce (has obsm["X_pca"],
        obsp["spatial_connectivities"], var["highly_variable"])
Output: AnnData with obs["spatial_cluster"], uns["moranI"],
        provenance in uns["omicsage_spatial_cluster"]

Two steps are combined here because SVG detection (Moran's I) requires
both the spatial neighbours graph (from spatial_reduce) and meaningful
cluster labels (from Leiden) for contextual interpretation.

Notes
-----
- Leiden clustering uses the KNN graph built from PCA (sc.pp.neighbors →
  sc.tl.leiden), NOT the spatial graph.  The spatial graph is used only
  for Moran's I.  This mirrors the sc-best-practices recommendation: use
  transcriptomic similarity for clustering, spatial context for SVG detection.
- Moran's I is run on HVGs only (to keep runtime tractable).
  Results are stored in uns["moranI"] (squidpy convention).
- annotation_map is optional: when provided, obs["spatial_cluster_label"]
  is written with human-readable names.  When None, only integer/string
  Leiden labels are written.
- Cluster numbering is non-deterministic — never hardcode cluster→celltype
  maps outside of a specific analysis run.
"""

from __future__ import annotations

from datetime import datetime
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spatial_cluster(
    adata: ad.AnnData,
    resolution: float = 0.5,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    random_state: int = 0,
    cluster_key: str = "spatial_cluster",
    annotation_map: Optional[dict] = None,
    run_svg: bool = True,
    svg_n_genes: Optional[int] = None,
    svg_n_jobs: int = 1,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Leiden clustering and Moran's I spatially variable gene detection.

    Parameters
    ----------
    adata
        AnnData produced by :func:`spatial_reduce`.  Must contain
        ``obsm["X_pca"]`` and ``obsp["spatial_connectivities"]``.
    resolution
        Leiden resolution parameter.  Higher values → more clusters.
        Default 0.5 is a reasonable starting point for Visium data.
    n_neighbors
        Number of transcriptomic neighbours for the KNN graph (used for
        Leiden clustering, not the spatial graph).
    n_pcs
        Number of PCA components to use when building the KNN graph.
    random_state
        Random seed for reproducibility.
    cluster_key
        Key written to ``obs`` for cluster labels.
        Default: ``"spatial_cluster"``.
    annotation_map
        Optional dict mapping Leiden cluster IDs (str) to human-readable
        cell type / region names.  When provided, writes
        ``obs[cluster_key + "_label"]``.  When ``None``, skipped.
    run_svg
        Whether to run Moran's I spatially variable gene detection.
        Requires ``obsp["spatial_connectivities"]`` to be present.
    svg_n_genes
        Number of HVGs to test for spatial autocorrelation.
        ``None`` (default) tests all HVGs.  Limit for speed on large datasets.
    svg_n_jobs
        Number of parallel jobs for Moran's I permutation test.
    inplace
        If ``False`` (default), operate on a copy.
        If ``True``, modify *adata* in place.

    Returns
    -------
    adata
        AnnData with:

        - ``obs[cluster_key]`` — Leiden cluster labels (str)
        - ``obs[cluster_key + "_label"]`` — human-readable labels
          (only if *annotation_map* provided)
        - ``uns["moranI"]`` — DataFrame of Moran's I results sorted by I
          (only if *run_svg* is ``True``)
        - ``uns["omicsage_spatial_cluster"]`` — provenance dict

    params
        Provenance dictionary.

    Raises
    ------
    ImportError
        If squidpy is not installed (required for Moran's I).
    TypeError
        If *adata* is not an AnnData object.
    ValueError
        If required keys are missing from *adata*.
    """
    _validate_input(adata, run_svg)

    adata = adata if inplace else adata.copy()

    # ------------------------------------------------------------------ #
    # 1.  KNN graph on PCA embedding (transcriptomic similarity)
    # ------------------------------------------------------------------ #
    n_pcs_actual = min(n_pcs, adata.obsm["X_pca"].shape[1])
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs_actual,
        random_state=random_state,
        use_rep="X_pca",
    )

    # ------------------------------------------------------------------ #
    # 2.  Leiden clustering
    # ------------------------------------------------------------------ #
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_state,
        key_added=cluster_key,
    )

    cluster_labels = adata.obs[cluster_key].astype(str)
    n_clusters = int(cluster_labels.nunique())
    cluster_sizes = cluster_labels.value_counts().to_dict()

    # ------------------------------------------------------------------ #
    # 3.  Optional annotation map → human-readable labels
    # ------------------------------------------------------------------ #
    label_key = cluster_key + "_label"
    if annotation_map is not None:
        adata.obs[label_key] = (
            cluster_labels.map(annotation_map).fillna("Unknown").astype("category")
        )
        n_annotated = int((adata.obs[label_key] != "Unknown").sum())
    else:
        n_annotated = 0

    # ------------------------------------------------------------------ #
    # 4.  Moran's I — spatially variable genes
    # ------------------------------------------------------------------ #
    svg_results: dict = {}
    n_svg_tested = 0

    if run_svg:
        _check_squidpy()
        # Select genes to test: HVGs only, optionally capped
        if "highly_variable" in adata.var.columns:
            hvg_mask = adata.var["highly_variable"].values
            gene_pool = adata.var_names[hvg_mask].tolist()
        else:
            gene_pool = adata.var_names.tolist()

        if svg_n_genes is not None:
            gene_pool = gene_pool[:svg_n_genes]

        n_svg_tested = len(gene_pool)

        sq.gr.spatial_autocorr(
            adata,
            mode="moran",
            genes=gene_pool,
            n_perms=None,   # analytical p-values (fast; no permutation)
            n_jobs=svg_n_jobs,
        )
        # Results written to adata.uns["moranI"] by squidpy
        if "moranI" in adata.uns:
            moran_df: pd.DataFrame = adata.uns["moranI"]
            n_sig = int((moran_df["pval_norm_fdr_bh"] < 0.05).sum())
            top5 = moran_df.head(5).index.tolist()
            svg_results = {
                "n_genes_tested": n_svg_tested,
                "n_significant_fdr05": n_sig,
                "top5_svg": top5,
            }

    # ------------------------------------------------------------------ #
    # 5.  Provenance
    # ------------------------------------------------------------------ #
    params = {
        "module": "spatial_cluster",
        "timestamp": datetime.now().isoformat(),
        "params": {
            "resolution": resolution,
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "n_pcs_actual": n_pcs_actual,
            "random_state": random_state,
            "cluster_key": cluster_key,
            "annotation_map_provided": annotation_map is not None,
            "run_svg": run_svg,
            "svg_n_genes": svg_n_genes,
        },
        "outputs": {
            "n_clusters": n_clusters,
            "cluster_sizes": cluster_sizes,
            "n_annotated_spots": n_annotated,
            **svg_results,
        },
    }
    adata.uns["omicsage_spatial_cluster"] = params
    return adata, params


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _check_squidpy() -> None:
    if not _SQUIDPY_AVAILABLE:
        raise ImportError(
            "squidpy is required for spatially variable gene detection. "
            "Install with: pip install squidpy"
        )


def _validate_input(adata: ad.AnnData, run_svg: bool) -> None:
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"Expected AnnData, got {type(adata).__name__}")
    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 observations (spots). Cannot cluster.")
    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 variables (genes). Cannot cluster.")
    if "X_pca" not in adata.obsm:
        raise ValueError(
            "adata.obsm['X_pca'] is missing. "
            "Run spatial_reduce() first."
        )
    if run_svg and "spatial_connectivities" not in adata.obsp:
        raise ValueError(
            "adata.obsp['spatial_connectivities'] is missing. "
            "Run spatial_reduce() first, or pass run_svg=False."
        )
