"""
spatial_qc.py — OmicSage Phase 7, Session 1
Visium spatial transcriptomics spot QC.

Input:  AnnData produced by spatial_ingest (has obsm["spatial"])
Output: AnnData with QC metrics, optional spot filtering, and
        provenance in uns["omicsage_spatial_qc"]
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


def spatial_qc(
    adata: ad.AnnData,
    min_counts: int = 500,
    max_counts: int = 100_000,
    min_genes: int = 200,
    max_genes: int = 10_000,
    max_mt_pct: float = 20.0,
    mt_prefix: str = "MT-",
    filter_spots: bool = True,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Compute QC metrics and optionally filter low-quality Visium spots.

    Parameters
    ----------
    adata
        AnnData produced by :func:`spatial_ingest`.  Must contain
        ``obsm["spatial"]``.
    min_counts
        Minimum total UMI counts per spot.
    max_counts
        Maximum total UMI counts per spot (removes ink/folding artefacts).
    min_genes
        Minimum number of detected genes per spot.
    max_genes
        Maximum number of detected genes per spot.
    max_mt_pct
        Maximum mitochondrial gene expression percentage per spot.
    mt_prefix
        Gene name prefix used to identify mitochondrial genes.
        ``"MT-"`` for human, ``"mt-"`` for mouse.
    filter_spots
        If ``True``, remove spots that fail any threshold.
        If ``False``, compute metrics only (thresholds stored in provenance).
    inplace
        If ``False`` (default), operate on a copy and return it.
        If ``True``, modify *adata* in place.

    Returns
    -------
    adata
        AnnData with QC columns in ``obs`` and provenance in
        ``uns["omicsage_spatial_qc"]``.
    params
        Dictionary of QC parameters, per-threshold counts, and summary stats.
    """
    _validate_input(adata)

    # Keep a reference to the caller's object for inplace provenance write.
    # filter_spots may rebind the local `adata` variable (new sliced copy),
    # which would break inplace=True if we don't track the original separately.
    original_adata = adata if inplace else None
    adata = adata if inplace else adata.copy()

    # ------------------------------------------------------------------ #
    # 1. Annotate mitochondrial genes
    # ------------------------------------------------------------------ #
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    n_mt_genes = int(adata.var["mt"].sum())

    # ------------------------------------------------------------------ #
    # 2. Calculate QC metrics (scanpy standard)
    # ------------------------------------------------------------------ #
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True,
        percent_top=None,
        log1p=False,
    )
    # Resulting obs columns:
    #   total_counts, n_genes_by_counts, pct_counts_mt

    # ------------------------------------------------------------------ #
    # 3. Compute per-threshold failure counts (before filtering)
    # ------------------------------------------------------------------ #
    n_spots_before = int(adata.n_obs)

    low_counts = int((adata.obs["total_counts"] < min_counts).sum())
    high_counts = int((adata.obs["total_counts"] > max_counts).sum())
    low_genes = int((adata.obs["n_genes_by_counts"] < min_genes).sum())
    high_genes = int((adata.obs["n_genes_by_counts"] > max_genes).sum())
    high_mt = int((adata.obs["pct_counts_mt"] > max_mt_pct).sum())

    # ------------------------------------------------------------------ #
    # 4. Build pass/fail mask
    # ------------------------------------------------------------------ #
    pass_mask = (
        (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["total_counts"] <= max_counts)
        & (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_mt_pct)
    )
    adata.obs["qc_pass"] = pass_mask

    n_spots_after = int(pass_mask.sum())
    n_spots_removed = n_spots_before - n_spots_after

    # ------------------------------------------------------------------ #
    # 5. Filter if requested
    # ------------------------------------------------------------------ #
    if filter_spots:
        adata = adata[pass_mask].copy()

    # ------------------------------------------------------------------ #
    # 6. Summary statistics (on retained spots)
    # ------------------------------------------------------------------ #
    summary_stats = _compute_summary_stats(adata)

    # ------------------------------------------------------------------ #
    # 7. Provenance
    # ------------------------------------------------------------------ #
    params = {
        "module": "spatial_qc",
        "timestamp": datetime.now().isoformat(),
        "params": {
            "min_counts": min_counts,
            "max_counts": max_counts,
            "min_genes": min_genes,
            "max_genes": max_genes,
            "max_mt_pct": max_mt_pct,
            "mt_prefix": mt_prefix,
            "filter_spots": filter_spots,
        },
        "outputs": {
            "n_spots_before": n_spots_before,
            "n_spots_after": n_spots_after,
            "n_spots_removed": n_spots_removed,
            "n_mt_genes_detected": n_mt_genes,
            "removed_low_counts": low_counts,
            "removed_high_counts": high_counts,
            "removed_low_genes": low_genes,
            "removed_high_genes": high_genes,
            "removed_high_mt": high_mt,
        },
        "summary_stats": summary_stats,
    }
    adata.uns["omicsage_spatial_qc"] = params
    # When inplace=True and filter_spots=True, `adata` is a new sliced copy
    # (adata[pass_mask].copy()), so we must also write provenance back to the
    # original input object the caller passed in.
    if inplace and original_adata is not None and original_adata is not adata:
        original_adata.uns["omicsage_spatial_qc"] = params
    return adata, params


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_input(adata: ad.AnnData) -> None:
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"Expected AnnData, got {type(adata).__name__}")
    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] is missing. "
            "Run spatial_ingest() first to load Visium data."
        )
    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 observations (spots). Cannot run QC.")
    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 variables (genes). Cannot run QC.")


def _compute_summary_stats(adata: ad.AnnData) -> dict:
    """Compute mean/median/std for key QC metrics on retained spots.

    Returns NaN for all stats when no spots remain after filtering.
    """
    stats = {}
    for col in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        if col in adata.obs.columns:
            vals = adata.obs[col].values
            if len(vals) == 0:
                stats[col] = {
                    "mean": float("nan"),
                    "median": float("nan"),
                    "std": float("nan"),
                    "min": float("nan"),
                    "max": float("nan"),
                }
            else:
                stats[col] = {
                    "mean": float(np.mean(vals)),
                    "median": float(np.median(vals)),
                    "std": float(np.std(vals)),
                    "min": float(np.min(vals)),
                    "max": float(np.max(vals)),
                }
    return stats
