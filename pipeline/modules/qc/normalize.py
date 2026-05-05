"""
OmicSage — Normalization Module
pipeline/modules/qc/normalize.py

Normalizes raw count AnnData (RNA) after QC filtering.
Input : AnnData with raw integer counts in .X  (output of qc.py → mdata["rna"])
Output: normalized AnnData  +  metrics dict

Steps
-----
1. Save raw counts to adata.layers['counts']
2. Normalize per-cell to target_sum counts  (scanpy normalize_total)
3. Log1p transform
4. Select top n_top_genes highly variable genes  (flavor='seurat_v3')
5. Store all parameters + software versions in adata.uns['omicsage_normalization']

Usage
-----
    from pipeline.modules.qc.normalize import normalize

    adata_rna = mdata["rna"]                          # pass RNA slot from run_qc()
    adata_norm, metrics = normalize(adata_rna)
    # or with custom params:
    adata_norm, metrics = normalize(adata_rna, target_sum=1e4, n_top_genes=3000)
"""

from __future__ import annotations

import logging
import warnings
from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def normalize(
    adata: AnnData,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    hvg_flavor: str = "seurat_v3",
    batch_key: Optional[str] = None,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
    inplace: bool = False,
    random_state: int = 0,
) -> tuple[AnnData, dict]:
    """
    Normalize raw-count RNA AnnData and select highly variable genes.

    Parameters
    ----------
    adata : AnnData
        AnnData with raw integer counts in .X.
        Typically ``mdata["rna"]`` from ``run_qc()``.
    target_sum : float
        Every cell is scaled to sum to this many counts before log1p.
        Default 1e4 (= "counts per 10 000", CP10K).
    n_top_genes : int
        Number of highly variable genes to flag.  Default 2000.
    hvg_flavor : str
        HVG selection method passed to ``sc.pp.highly_variable_genes``.
        'seurat_v3' (default) uses the variance-stabilization approach on
        raw counts, so it is run *before* log1p transform.
        Use 'seurat' or 'cell_ranger' only if counts are already in .X after log1p.
    batch_key : str, optional
        Column in ``adata.obs`` to use as batch label during HVG selection.
        When provided, HVGs are selected per batch and the results are
        combined (genes highly variable in at least one batch are flagged).
        Recommended for multi-donor / multi-site data — prevents a single
        dominant batch from monopolising the HVG list.
        Typical values: ``'batch'``, ``'Site'``, ``'DonorID'``.
        Default None (HVG selection ignores batch structure).
    min_mean, max_mean, min_disp : float
        Additional HVG filters (only active for 'seurat' / 'cell_ranger' flavors).
    inplace : bool
        If False (default) work on a copy so the caller's object is not mutated.
    random_state : int
        Passed to scanpy functions for reproducibility.

    Returns
    -------
    adata_norm : AnnData
        Normalized, log-transformed AnnData.
        - ``.X``                          — log1p-normalized values
        - ``.layers['counts']``           — original raw integer counts
        - ``.var['highly_variable']``     — boolean HVG flag
        - ``.uns['omicsage_normalization']`` — full parameter + metrics record
    metrics : dict
        Summary statistics for logging / downstream QA.
    """
    _validate_input(adata)

    # Work on a copy unless caller opts in to inplace mutation
    adata_out = adata if inplace else adata.copy()

    # ------------------------------------------------------------------
    # Step 1 — Persist raw counts before any modification
    # ------------------------------------------------------------------
    if "counts" not in adata_out.layers:
        logger.info("Saving raw counts to adata.layers['counts']")
        adata_out.layers["counts"] = adata_out.X.copy()
    else:
        logger.info("adata.layers['counts'] already exists — skipping copy")

    # ------------------------------------------------------------------
    # Step 2 — HVG selection on raw counts (seurat_v3 requirement)
    # seurat_v3 expects counts, not normalized values, so we do it first.
    # For other flavors HVG will be re-run after normalization below.
    # ------------------------------------------------------------------
    if hvg_flavor == "seurat_v3":
        logger.info(
            "Selecting top %d HVGs with flavor='seurat_v3' on raw counts%s",
            n_top_genes,
            f" (batch_key='{batch_key}')" if batch_key else "",
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.highly_variable_genes(
                adata_out,
                flavor=hvg_flavor,
                n_top_genes=n_top_genes,
                layer="counts",      # explicitly use raw counts layer
                batch_key=batch_key,
            )

    # ------------------------------------------------------------------
    # Step 3 — Normalize per cell to target_sum
    # ------------------------------------------------------------------
    logger.info("Normalizing to target_sum=%.0f", target_sum)
    sc.pp.normalize_total(adata_out, target_sum=target_sum)

    # ------------------------------------------------------------------
    # Step 4 — Log1p transform
    # ------------------------------------------------------------------
    logger.info("Applying log1p transform")
    sc.pp.log1p(adata_out)

    # Save log1p-normalized values to layers so both raw and normalized
    # counts are always accessible without recomputation
    adata_out.layers["logcounts"] = adata_out.X.copy()
    logger.info("Saved log1p-normalized values to adata.layers['logcounts']")

    # ------------------------------------------------------------------
    # Step 5 — HVG for non-seurat_v3 flavors (run post log1p)
    # ------------------------------------------------------------------
    if hvg_flavor != "seurat_v3":
        logger.info(
            "Selecting top %d HVGs with flavor='%s' on log-normalized data%s",
            n_top_genes,
            hvg_flavor,
            f" (batch_key='{batch_key}')" if batch_key else "",
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.highly_variable_genes(
                adata_out,
                flavor=hvg_flavor,
                n_top_genes=n_top_genes,
                min_mean=min_mean,
                max_mean=max_mean,
                min_disp=min_disp,
                batch_key=batch_key,
            )

    # ------------------------------------------------------------------
    # Step 6 — Build metrics dict
    # ------------------------------------------------------------------
    n_hvg_selected = int(adata_out.var["highly_variable"].sum())
    x_array = (
        adata_out.X.toarray() if hasattr(adata_out.X, "toarray") else np.asarray(adata_out.X)
    )
    metrics = {
        "n_cells": int(adata_out.n_obs),
        "n_genes": int(adata_out.n_vars),
        "n_hvg_selected": n_hvg_selected,
        "target_sum": target_sum,
        "hvg_flavor": hvg_flavor,
        "batch_key": batch_key,
        "mean_counts_per_cell_after_norm": float(x_array.sum(axis=1).mean()),
        "log1p_applied": True,
        "raw_counts_in_layer": "counts",
        "normalized_in_layer": "logcounts",
    }

    # ------------------------------------------------------------------
    # Step 7 — Store provenance in uns
    # ------------------------------------------------------------------
    import scanpy
    adata_out.uns["omicsage_normalization"] = {
        "target_sum": target_sum,
        "n_top_genes": n_top_genes,
        "hvg_flavor": hvg_flavor,
        "batch_key": batch_key,
        "min_mean": min_mean,
        "max_mean": max_mean,
        "min_disp": min_disp,
        "n_hvg_selected": n_hvg_selected,
        "log1p_applied": True,
        "normalized_in_layer": "logcounts",
        "scanpy_version": scanpy.__version__,
        "omicsage_module": "pipeline.modules.qc.normalize",
        "omicsage_version": "0.1.0",
    }

    logger.info(
        "Normalization complete — %d cells, %d genes, %d HVGs selected",
        adata_out.n_obs,
        adata_out.n_vars,
        n_hvg_selected,
    )

    return adata_out, metrics


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_input(adata: AnnData) -> None:
    """Raise informative errors for common input mistakes."""
    if not isinstance(adata, AnnData):
        raise TypeError(
            f"normalize() expects an AnnData object, got {type(adata).__name__}. "
            "Pass mdata['rna'] not the full MuData."
        )

    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 cells — nothing to normalize.")

    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 genes — nothing to normalize.")

    # Check X looks like raw counts (integers or very close to integers)
    x_sample = adata.X[:100, :100]
    if hasattr(x_sample, "toarray"):
        x_sample = x_sample.toarray()
    x_sample = np.asarray(x_sample, dtype=float)

    non_integer_frac = np.mean(np.abs(x_sample - np.round(x_sample)) > 1e-3)
    if non_integer_frac > 0.05:
        raise ValueError(
            "More than 5% of sampled values in adata.X are non-integer. "
            "normalize() expects raw counts in .X. "
            "If you already normalized, do not call normalize() again. "
            f"(Non-integer fraction sampled: {non_integer_frac:.2%})"
        )

    logger.debug("Input validation passed — %d cells × %d genes", adata.n_obs, adata.n_vars)
