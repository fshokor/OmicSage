"""
spatial_downstream.py — OmicSage Phase 7, Session 5
Spatial downstream analyses: region clustering, cell-type specific expression,
cell-type specific SVGs, co-occurrence, neighbourhood enrichment,
ligand-receptor communication, and pathway enrichment on SVGs.

Input:  AnnData from spatial_deconvolve (has obsm["q05_cell_abundance_w_sf"],
        obs["dominant_cell_type"], uns["moranI"])
Output: AnnData with downstream analysis results in uns, provenance in
        uns["omicsage_spatial_downstream"]

All analyses gracefully skip if their required inputs are absent.

Output contract
---------------
uns["region_cluster_results"]        — region clustering provenance
uns["celltype_marker_genes"]         — dict: cell_type → list of top correlated genes
uns["celltype_svg"]                  — dict: cell_type → DataFrame (Moran's I per subset)
uns["{dominant_celltype_key}_co_occurrence"]      — set by sq.gr.co_occurrence
uns["{dominant_celltype_key}_nhood_enrichment"]   — set by sq.gr.nhood_enrichment
uns["{dominant_celltype_key}_ligrec"]             — set by sq.gr.ligrec
uns["svg_gsea"]                      — DataFrame: top enriched pathways
uns["omicsage_spatial_downstream"]   — provenance dict

obs["region_cluster"]                — Leiden clusters from cell type abundance
obsm["X_umap_celltype"]             — UMAP of region clusters
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

logger = logging.getLogger(__name__)

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False

try:
    import gseapy
    _GSEAPY_AVAILABLE = True
except ImportError:
    _GSEAPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spatial_downstream(
    adata: ad.AnnData,
    # Region clustering from cell type abundance (sc-best-practices ch.32.3.4.2)
    run_region_clustering: bool = True,
    region_resolution: float = 0.5,
    region_n_neighbors: int = 15,
    # Cell-type specific gene expression (Spearman)
    run_celltype_expression: bool = True,
    n_marker_genes: int = 20,
    # Cell-type specific SVGs
    run_celltype_svg: bool = True,
    svg_n_genes: Optional[int] = None,
    # Co-occurrence
    run_co_occurrence: bool = True,
    co_occurrence_interval: Optional[list] = None,
    # Neighbourhood enrichment
    run_nhood_enrichment: bool = True,
    n_perms_nhood: int = 1000,
    # Ligand-receptor
    run_ligrec: bool = True,
    ligrec_n_perms: int = 1000,
    ligrec_organism: str = "human",
    # SVG pathway enrichment
    run_svg_gsea: bool = True,
    svg_gsea_gene_sets: str = "GO_Biological_Process_2023",
    svg_gsea_organism: str = "Human",
    # Cell type column (from deconvolution)
    dominant_celltype_key: str = "dominant_cell_type",
    # General
    n_jobs: int = 1,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Run all spatial downstream analyses.

    Parameters
    ----------
    adata
        AnnData produced by :func:`spatial_deconvolve`. Should contain
        ``obsm["q05_cell_abundance_w_sf"]``, ``obs["dominant_cell_type"]``,
        ``uns["moranI"]``, and ``obsp["spatial_connectivities"]``.
        Each analysis gracefully skips if its required keys are absent.
    run_region_clustering
        Cluster spots by cell type composition (sc-best-practices ch.32.3.4.2).
        Requires ``obsm["q05_cell_abundance_w_sf"]``.
    region_resolution
        Leiden resolution for region clustering. Default 0.5.
    region_n_neighbors
        Neighbours for region KNN graph. Capped at n_cell_types - 1.
    run_celltype_expression
        Spearman-correlate each cell type's abundance against all genes.
        Requires ``obsm["q05_cell_abundance_w_sf"]``.
    n_marker_genes
        Top N correlated genes to store per cell type.
    run_celltype_svg
        Moran's I on spots enriched for each cell type.
        Requires ``uns["moranI"]`` + ``obsm["q05_cell_abundance_w_sf"]``.
    svg_n_genes
        Max genes to test per cell type subset. ``None`` = all HVGs.
    run_co_occurrence
        ``sq.gr.co_occurrence`` — co-localisation across distance ranges.
        Requires ``obs[dominant_celltype_key]``.
    co_occurrence_interval
        Interval range list. ``None`` = squidpy default (50 intervals).
    run_nhood_enrichment
        ``sq.gr.nhood_enrichment`` — permutation test on spatial adjacency.
        Requires ``obsp["spatial_connectivities"]`` + ``obs[dominant_celltype_key]``.
    n_perms_nhood
        Permutations for neighbourhood enrichment. Default 1000.
    run_ligrec
        ``sq.gr.ligrec`` — ligand-receptor communication. Requires gene symbols
        in ``var["feature_name"]`` or as ``var_names``.
    ligrec_n_perms
        Permutations for ligrec. Default 1000.
    ligrec_organism
        Organism for OmniPath database lookup: ``"human"`` or ``"mouse"``.
    run_svg_gsea
        ``gseapy.prerank`` on the Moran's I ranked gene list.
        Requires ``uns["moranI"]`` + gseapy installed.
    svg_gsea_gene_sets
        Gene set collection name (Enrichr) or path to .gmt file.
    svg_gsea_organism
        Organism for Enrichr lookup: ``"Human"`` or ``"Mouse"``.
    dominant_celltype_key
        ``obs`` column with dominant cell type per spot.
    n_jobs
        Parallel workers for analyses that support it.
    inplace
        If ``False`` (default), operate on a copy.

    Returns
    -------
    adata
        AnnData with downstream results written to ``uns``.
    params
        Provenance dictionary (same as ``uns["omicsage_spatial_downstream"]``).
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"Expected AnnData, got {type(adata).__name__}")

    adata = adata if inplace else adata.copy()

    provenance: dict = {
        "module": "spatial_downstream",
        "timestamp": datetime.now().isoformat(),
        "params": {
            "run_region_clustering": run_region_clustering,
            "region_resolution": region_resolution,
            "run_celltype_expression": run_celltype_expression,
            "n_marker_genes": n_marker_genes,
            "run_celltype_svg": run_celltype_svg,
            "run_co_occurrence": run_co_occurrence,
            "run_nhood_enrichment": run_nhood_enrichment,
            "n_perms_nhood": n_perms_nhood,
            "run_ligrec": run_ligrec,
            "ligrec_organism": ligrec_organism,
            "run_svg_gsea": run_svg_gsea,
            "svg_gsea_gene_sets": svg_gsea_gene_sets,
            "dominant_celltype_key": dominant_celltype_key,
        },
        "analyses": {},
    }

    cell_types = _get_cell_types(adata, dominant_celltype_key)

    # ------------------------------------------------------------------ #
    # 0. Region clustering from cell type abundance
    # ------------------------------------------------------------------ #
    if run_region_clustering:
        print("  [downstream] 0/6 region clustering ...", flush=True)
        result = _run_region_clustering(
            adata, region_resolution, region_n_neighbors
        )
        provenance["analyses"]["region_clustering"] = result
        if not result.get("skipped"):
            logger.info(
                "[downstream] region clustering: %d regions",
                result.get("n_regions", 0),
            )

    # ------------------------------------------------------------------ #
    # 1. Cell-type specific gene expression (Spearman correlation)
    # ------------------------------------------------------------------ #
    if run_celltype_expression:
        print("  [downstream] 1/6 cell-type expression ...", flush=True)
        result = _run_celltype_expression(
            adata, cell_types, n_marker_genes
        )
        provenance["analyses"]["celltype_expression"] = result
        if not result.get("skipped"):
            logger.info(
                "[downstream] celltype expression: %d cell types profiled",
                result.get("n_cell_types", 0),
            )

    # ------------------------------------------------------------------ #
    # 2. Cell-type specific SVGs
    # ------------------------------------------------------------------ #
    if run_celltype_svg:
        print("  [downstream] 2/6 cell-type SVGs ...", flush=True)
        result = _run_celltype_svg(
            adata, cell_types, svg_n_genes, n_jobs
        )
        provenance["analyses"]["celltype_svg"] = result
        if not result.get("skipped"):
            logger.info(
                "[downstream] celltype SVGs: %d cell types",
                result.get("n_cell_types", 0),
            )

    # ------------------------------------------------------------------ #
    # 3. Co-occurrence
    # ------------------------------------------------------------------ #
    if run_co_occurrence:
        print("  [downstream] 3/6 co-occurrence ...", flush=True)
        result = _run_co_occurrence(
            adata, dominant_celltype_key, co_occurrence_interval, n_jobs
        )
        provenance["analyses"]["co_occurrence"] = result
        if not result.get("skipped"):
            logger.info("[downstream] co-occurrence complete")

    # ------------------------------------------------------------------ #
    # 4. Neighbourhood enrichment
    # ------------------------------------------------------------------ #
    if run_nhood_enrichment:
        print("  [downstream] 4/6 neighbourhood enrichment ...", flush=True)
        result = _run_nhood_enrichment(
            adata, dominant_celltype_key, n_perms_nhood, n_jobs
        )
        provenance["analyses"]["nhood_enrichment"] = result
        if not result.get("skipped"):
            logger.info("[downstream] nhood enrichment complete")

    # ------------------------------------------------------------------ #
    # 5. Ligand-receptor communication
    # ------------------------------------------------------------------ #
    if run_ligrec:
        print("  [downstream] 5/6 ligand-receptor ...", flush=True)
        result = _run_ligrec(
            adata, dominant_celltype_key, ligrec_n_perms, ligrec_organism, n_jobs
        )
        provenance["analyses"]["ligrec"] = result
        if not result.get("skipped"):
            logger.info("[downstream] ligrec complete")

    # ------------------------------------------------------------------ #
    # 6. SVG pathway enrichment
    # ------------------------------------------------------------------ #
    if run_svg_gsea:
        print("  [downstream] 6/6 SVG pathway enrichment ...", flush=True)
        result = _run_svg_gsea(
            adata, svg_gsea_gene_sets, svg_gsea_organism, n_jobs
        )
        provenance["analyses"]["svg_gsea"] = result
        if not result.get("skipped"):
            logger.info(
                "[downstream] SVG GSEA: %d pathways, %d significant",
                result.get("n_pathways", 0),
                result.get("n_significant", 0),
            )

    adata.uns["omicsage_spatial_downstream"] = provenance
    return adata, provenance


# ---------------------------------------------------------------------------
# Internal helpers — cell type resolution
# ---------------------------------------------------------------------------


def _get_cell_types(adata: ad.AnnData, dominant_celltype_key: str) -> list[str]:
    """Return list of cell type names from deconvolution provenance.

    cell_type_names may be deserialized from h5ad as a numpy array — coerce to
    list before any truth-value check to avoid the ambiguous array bool error.
    """
    deconv_prov = adata.uns.get("omicsage_spatial_deconvolve", {})
    cell_types = list(deconv_prov.get("outputs", {}).get("cell_type_names", []))
    if len(cell_types) > 0:
        return [ct for ct in cell_types if ct in adata.obs.columns]

    # Fallback: numeric obs columns that are not standard QC / pipeline keys
    _skip = {
        "total_counts", "n_genes_by_counts", "log1p_total_counts",
        "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_50_genes",
        "qc_pass", dominant_celltype_key, "region_cluster", "leiden",
        "spatial_cluster", "spatial_cluster_label",
        "array_row", "array_col", "in_tissue",
    }
    return [
        c for c in adata.obs.columns
        if c not in _skip and pd.api.types.is_numeric_dtype(adata.obs[c])
    ]


# ---------------------------------------------------------------------------
# Analysis 0 — Region clustering from cell type abundance
# ---------------------------------------------------------------------------


def _run_region_clustering(
    adata: ad.AnnData,
    resolution: float,
    n_neighbors: int,
) -> dict:
    """Cluster spots by cell type composition (sc-best-practices §32.3.4.2)."""
    if "q05_cell_abundance_w_sf" not in adata.obsm:
        return {"skipped": True, "reason": "obsm['q05_cell_abundance_w_sf'] not found"}

    try:
        abundance = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])
        n_dims = abundance.shape[1]
        k = min(n_neighbors, n_dims - 1, adata.n_obs - 1)

        sc.pp.neighbors(
            adata,
            use_rep="q05_cell_abundance_w_sf",
            n_neighbors=k,
            key_added="neighbors_celltype",
        )
        sc.tl.leiden(
            adata,
            neighbors_key="neighbors_celltype",
            key_added="region_cluster",
            resolution=resolution,
        )
        adata.obs["region_cluster"] = adata.obs["region_cluster"].astype("category")
        n_regions = int(adata.obs["region_cluster"].nunique())

        # UMAP for visualization — save under non-default key to avoid overwriting
        # gene-expression UMAP produced by other pipeline steps
        try:
            _prev_umap = adata.obsm.pop("X_umap", None)
            sc.tl.umap(adata, neighbors_key="neighbors_celltype", min_dist=0.3)
            adata.obsm["X_umap_celltype"] = adata.obsm.pop("X_umap")
            if _prev_umap is not None:
                adata.obsm["X_umap"] = _prev_umap
        except Exception as e_umap:
            logger.warning("[downstream] UMAP for region clustering failed: %s", e_umap)

        return {"skipped": False, "n_regions": n_regions, "resolution": resolution}

    except Exception as e:
        return {"skipped": True, "reason": f"region clustering failed: {e}"}


# ---------------------------------------------------------------------------
# Analysis 1 — Cell-type specific gene expression
# ---------------------------------------------------------------------------


def _run_celltype_expression(
    adata: ad.AnnData,
    cell_types: list[str],
    n_marker_genes: int,
) -> dict:
    """Spearman correlation between cell type abundance and gene expression."""
    if "q05_cell_abundance_w_sf" not in adata.obsm:
        return {"skipped": True, "reason": "obsm['q05_cell_abundance_w_sf'] not found"}
    if not cell_types:
        return {"skipped": True, "reason": "no cell types found"}

    try:
        from scipy.stats import spearmanr

        # Expression matrix — work column-by-column to avoid densifying
        X = adata.X
        is_sparse = sp.issparse(X)

        marker_dict: dict[str, list[str]] = {}

        for ct in cell_types:
            if ct not in adata.obs.columns:
                continue
            abundance = adata.obs[ct].values.astype(np.float32)

            # Compute Spearman r for every gene using dense columns
            # (safe for typical Visium: ~33k genes, ~10k spots)
            if is_sparse:
                expr_dense = np.asarray(X.toarray(), dtype=np.float32)
            else:
                expr_dense = np.asarray(X, dtype=np.float32)

            # Vectorised rank-based Spearman (avoids per-gene scipy calls)
            from scipy.stats import rankdata
            ranked_abundance = rankdata(abundance).astype(np.float32)
            ra = ranked_abundance - ranked_abundance.mean()

            ranked_expr = np.apply_along_axis(rankdata, 0, expr_dense).astype(np.float32)
            re = ranked_expr - ranked_expr.mean(axis=0)

            numer = (ra[:, None] * re).sum(axis=0)
            denom = (
                np.sqrt((ra ** 2).sum())
                * np.sqrt((re ** 2).sum(axis=0))
            )
            denom = np.where(denom == 0, 1e-10, denom)
            corr = numer / denom

            top_idx = np.argsort(np.abs(corr))[::-1][:n_marker_genes]
            top_genes = adata.var_names[top_idx].tolist()

            # Prefer gene symbols when available
            if "feature_name" in adata.var.columns:
                symbol_map = adata.var["feature_name"].to_dict()
                top_genes = [symbol_map.get(g, g) for g in top_genes]

            marker_dict[ct] = top_genes

        adata.uns["celltype_marker_genes"] = marker_dict

        return {
            "skipped": False,
            "n_cell_types": len(marker_dict),
            "n_marker_genes": n_marker_genes,
        }

    except Exception as e:
        return {"skipped": True, "reason": f"celltype expression failed: {e}"}


# ---------------------------------------------------------------------------
# Analysis 2 — Cell-type specific SVGs
# ---------------------------------------------------------------------------


def _run_celltype_svg(
    adata: ad.AnnData,
    cell_types: list[str],
    svg_n_genes: Optional[int],
    n_jobs: int,
) -> dict:
    """Moran's I on spots enriched for each cell type."""
    if not _SQUIDPY_AVAILABLE:
        return {"skipped": True, "reason": "squidpy not installed"}
    if "moranI" not in adata.uns:
        return {"skipped": True, "reason": "uns['moranI'] not found — run spatial_cluster first"}
    if not cell_types:
        return {"skipped": True, "reason": "no cell types found"}
    if "q05_cell_abundance_w_sf" not in adata.obsm:
        return {"skipped": True, "reason": "obsm['q05_cell_abundance_w_sf'] not found"}

    celltype_svg: dict[str, pd.DataFrame] = {}

    for ct in cell_types:
        if ct not in adata.obs.columns:
            continue
        ct_median = float(adata.obs[ct].median())
        mask = adata.obs[ct].values > ct_median
        if mask.sum() < 10:
            logger.debug("[downstream] skipping cell type SVG for %s: < 10 spots", ct)
            continue

        try:
            ad_sub = adata[mask].copy()

            # Recompute spatial neighbours on the subset (original graph is invalid)
            sq.gr.spatial_neighbors(
                ad_sub, n_neighs=6, coord_type=None, key_added="spatial"
            )

            # Select genes
            if "highly_variable" in ad_sub.var.columns:
                gene_pool = ad_sub.var_names[
                    ad_sub.var["highly_variable"].values
                ].tolist()
            else:
                gene_pool = ad_sub.var_names.tolist()

            if svg_n_genes is not None:
                gene_pool = gene_pool[:svg_n_genes]

            if not gene_pool:
                continue

            sq.gr.spatial_autocorr(
                ad_sub,
                mode="moran",
                genes=gene_pool,
                n_perms=None,   # analytical p-values for speed
                n_jobs=n_jobs,
            )

            if "moranI" in ad_sub.uns:
                celltype_svg[ct] = ad_sub.uns["moranI"]

        except Exception as e:
            logger.warning("[downstream] celltype SVG failed for %s: %s", ct, e)

    adata.uns["celltype_svg"] = celltype_svg

    return {
        "skipped": False,
        "n_cell_types": len(celltype_svg),
    }


# ---------------------------------------------------------------------------
# Analysis 3 — Co-occurrence
# ---------------------------------------------------------------------------


def _run_co_occurrence(
    adata: ad.AnnData,
    dominant_celltype_key: str,
    interval: Optional[list],
    n_jobs: int,
) -> dict:
    """sq.gr.co_occurrence across spatial distance ranges."""
    if not _SQUIDPY_AVAILABLE:
        return {"skipped": True, "reason": "squidpy not installed"}
    if dominant_celltype_key not in adata.obs.columns:
        return {
            "skipped": True,
            "reason": f"obs['{dominant_celltype_key}'] not found",
        }
    if "spatial" not in adata.obsm:
        return {"skipped": True, "reason": "obsm['spatial'] not found"}

    try:
        kwargs: dict = {
            "cluster_key": dominant_celltype_key,
            "n_jobs": n_jobs,
        }
        if interval is not None:
            kwargs["interval"] = interval

        sq.gr.co_occurrence(adata, **kwargs)

        return {"skipped": False}

    except Exception as e:
        return {"skipped": True, "reason": f"co-occurrence failed: {e}"}


# ---------------------------------------------------------------------------
# Analysis 4 — Neighbourhood enrichment
# ---------------------------------------------------------------------------


def _run_nhood_enrichment(
    adata: ad.AnnData,
    dominant_celltype_key: str,
    n_perms: int,
    n_jobs: int,
) -> dict:
    """sq.gr.nhood_enrichment — permutation test on spatial adjacency."""
    if not _SQUIDPY_AVAILABLE:
        return {"skipped": True, "reason": "squidpy not installed"}
    if dominant_celltype_key not in adata.obs.columns:
        return {
            "skipped": True,
            "reason": f"obs['{dominant_celltype_key}'] not found",
        }
    if "spatial_connectivities" not in adata.obsp:
        # Try to recompute spatial graph as a safety net
        if "spatial" in adata.obsm:
            logger.info(
                "[downstream] spatial_connectivities missing; recomputing graph"
            )
            try:
                sq.gr.spatial_neighbors(
                    adata, n_neighs=6, coord_type=None, key_added="spatial"
                )
            except Exception as e:
                return {
                    "skipped": True,
                    "reason": f"spatial graph recompute failed: {e}",
                }
        else:
            return {
                "skipped": True,
                "reason": "obsp['spatial_connectivities'] not found",
            }

    try:
        sq.gr.nhood_enrichment(
            adata,
            cluster_key=dominant_celltype_key,
            n_perms=n_perms,
            n_jobs=n_jobs,
            seed=42,
        )
        return {"skipped": False, "n_perms": n_perms}

    except Exception as e:
        return {"skipped": True, "reason": f"nhood enrichment failed: {e}"}


# ---------------------------------------------------------------------------
# Analysis 5 — Ligand-receptor communication
# ---------------------------------------------------------------------------


def _serialize_ligrec_uns(adata: ad.AnnData, ligrec_key: str) -> None:
    """Serialize squidpy ligrec uns entry to h5ad-compatible format.

    sq.gr.ligrec stores three DataFrames under uns[ligrec_key]:
        means     — mean expression per LR pair per cell-type pair
        pvalues   — permutation p-values (same shape)
        metadata  — LR pair metadata (gene names, categories)

    Each DataFrame has a pd.MultiIndex on both axes (LR pairs × cell-type pairs)
    which anndata cannot write to HDF5. We serialize each DataFrame as a dict:
        {
          "<key>_data":    df.to_json(orient="records"),   # row data
          "<key>_columns": json.dumps([str(c) for c in df.columns]),
          "<key>_index":   json.dumps([str(i) for i in df.index]),
        }
    This round-trips perfectly regardless of shape.
    """
    import json as _json
    raw = adata.uns[ligrec_key]
    serialized: dict = {}
    for key, val in raw.items():
        if isinstance(val, pd.DataFrame):
            df = val.copy()
            # Stringify both axes before serializing
            df.index   = [str(i) for i in df.index]
            df.columns = [str(c) for c in df.columns]
            serialized[f"{key}_data"]    = df.to_json(orient="records")
            serialized[f"{key}_columns"] = _json.dumps(list(df.columns))
            serialized[f"{key}_index"]   = _json.dumps(list(df.index))
        else:
            try:
                serialized[key] = val
            except Exception:
                serialized[key] = str(val)
    adata.uns[ligrec_key] = serialized


def _run_ligrec(
    adata: ad.AnnData,
    dominant_celltype_key: str,
    n_perms: int,
    organism: str,
    n_jobs: int,
) -> dict:
    """sq.gr.ligrec with gene symbol swap for ENSEMBL-indexed data."""
    if not _SQUIDPY_AVAILABLE:
        return {"skipped": True, "reason": "squidpy not installed"}
    if not hasattr(sq.gr, "ligrec"):
        return {"skipped": True, "reason": "sq.gr.ligrec not available in this squidpy version"}
    if dominant_celltype_key not in adata.obs.columns:
        return {
            "skipped": True,
            "reason": f"obs['{dominant_celltype_key}'] not found",
        }

    try:
        # ligrec matches gene names against OmniPath database (gene symbols).
        # When var_names are ENSEMBL IDs, we must temporarily remap to symbols.
        if "feature_name" in adata.var.columns:
            ad_sym = adata.copy()
            ad_sym.var_names = (
                ad_sym.var["feature_name"].astype(str).values
            )
            ad_sym.var_names_make_unique()
            sq.gr.ligrec(
                ad_sym,
                n_perms=n_perms,
                cluster_key=dominant_celltype_key,
                use_raw=False,
                n_jobs=n_jobs,
                seed=42,
            )
            ligrec_key = f"{dominant_celltype_key}_ligrec"
            if ligrec_key in ad_sym.uns:
                adata.uns[ligrec_key] = ad_sym.uns[ligrec_key]
            del ad_sym
        else:
            sq.gr.ligrec(
                adata,
                n_perms=n_perms,
                cluster_key=dominant_celltype_key,
                use_raw=False,
                n_jobs=n_jobs,
                seed=42,
            )

        ligrec_key = f"{dominant_celltype_key}_ligrec"
        has_result = ligrec_key in adata.uns
        if has_result:
            _serialize_ligrec_uns(adata, ligrec_key)
        return {"skipped": not has_result, "reason": "no result stored" if not has_result else None}

    except Exception as e:
        return {"skipped": True, "reason": f"ligrec failed: {e}"}


# ---------------------------------------------------------------------------
# Analysis 6 — SVG pathway enrichment
# ---------------------------------------------------------------------------

# gseapy.prerank returns these columns as object dtype in most versions;
# they must be cast to float64 before anndata can serialise them to h5ad.
_GSEA_FLOAT_COLS = {"ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"}


def _sanitize_gsea_df(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce gseapy res2d to h5ad-serialisable dtypes.

    Numeric columns (ES, NES, p-values) are cast to float64.
    All other columns (Term, Tag %, Lead_genes, …) are cast to str.
    This prevents the ``TypeError: Can't implicitly convert non-string
    objects to strings`` crash when anndata writes object-dtype columns
    to HDF5.
    """
    df = df.copy()
    for col in df.columns:
        if col in _GSEA_FLOAT_COLS:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")
        else:
            df[col] = df[col].astype(str)
    return df


def _run_svg_gsea(
    adata: ad.AnnData,
    gene_sets: str,
    organism: str,
    n_jobs: int,
) -> dict:
    """gseapy.prerank on Moran's I ranked SVG list."""
    if not _GSEAPY_AVAILABLE:
        return {"skipped": True, "reason": "gseapy not installed"}
    if "moranI" not in adata.uns:
        return {"skipped": True, "reason": "uns['moranI'] not found — run spatial_cluster first"}

    try:
        import gseapy

        moranI_df = adata.uns["moranI"].copy()

        # Map ENSEMBL IDs → gene symbols when available
        if "feature_name" in adata.var.columns:
            symbol_map = adata.var["feature_name"].to_dict()
            moranI_df.index = [symbol_map.get(g, g) for g in moranI_df.index]

        ranked = moranI_df["I"].sort_values(ascending=False)
        # Clean: remove missing names and duplicates
        ranked = ranked[ranked.index.notna()]
        ranked = ranked[ranked.index != ""]
        ranked = ranked[~ranked.index.duplicated(keep="first")]

        enr = gseapy.prerank(
            rnk=ranked,
            gene_sets=gene_sets,
            outdir=None,
            no_plot=True,
            min_size=5,
            max_size=500,
            permutation_num=100,
            threads=n_jobs,
            organism=organism,
            verbose=False,
        )

        # Sanitise before storing: gseapy returns numeric cols as object dtype
        # which anndata cannot write to HDF5 without explicit coercion.
        gsea_df = _sanitize_gsea_df(enr.res2d)
        adata.uns["svg_gsea"] = gsea_df

        n_pathways = len(gsea_df)
        fdr_col = next(
            (c for c in ["FDR q-val", "fdr"] if c in gsea_df.columns), None
        )
        n_sig = int((gsea_df[fdr_col] < 0.05).sum()) if fdr_col else 0

        return {
            "skipped": False,
            "n_pathways": n_pathways,
            "n_significant": n_sig,
            "gene_sets": gene_sets,
        }

    except Exception as e:
        return {"skipped": True, "reason": f"SVG GSEA failed: {e}"}
