"""
OmicSage — QC Module
pipeline/modules/qc/qc.py

Input:  AnnData with raw counts in adata.X (output of ingest.py)
Output: Filtered AnnData + QC metrics dict

Metrics computed per cell:
  - n_genes_by_counts   : number of genes detected
  - total_counts        : total UMI count
  - pct_counts_mt       : mitochondrial read percentage

Filters applied (all configurable):
  - min_genes           : minimum genes per cell (default 200)
  - max_genes           : maximum genes per cell (default 6000)
  - max_mt_pct          : maximum MT% per cell   (default 20)
  - remove_doublets     : drop Scrublet doublets  (default True)
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import scanpy as sc
import scrublet as scr
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_qc(
    adata: AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_mt_pct: float = 20.0,
    remove_doublets: bool = True,
    scrublet_expected_doublet_rate: float = 0.06,
    random_state: int = 0,
    generate_report: bool = False,
    report_path: str = "reports/qc_report.html",
    sample_name: str = "sample",
) -> tuple[AnnData, dict[str, Any]]:
    """Run QC on a raw-count AnnData object.

    Parameters
    ----------
    adata:
        AnnData with raw integer counts in ``adata.X``.
        Produced by ``ingest.py``.
    min_genes:
        Minimum number of genes detected per cell.
    max_genes:
        Maximum number of genes detected per cell
        (removes likely doublets / multiplets not caught by Scrublet).
    max_mt_pct:
        Maximum mitochondrial read percentage per cell.
    remove_doublets:
        If True, remove cells flagged as doublets by Scrublet.
    scrublet_expected_doublet_rate:
        Prior doublet rate passed to Scrublet (default 6 %).
    random_state:
        Reproducibility seed for Scrublet.
    generate_report:
        If True, write a self-contained HTML QC report after filtering.
    report_path:
        Output path for the HTML report (used only if generate_report=True).
    sample_name:
        Label shown in the report header.

    Returns
    -------
    adata_filtered:
        AnnData containing only cells that passed all QC filters.
        QC columns are added to ``adata_filtered.obs``.
    metrics:
        Summary dict with counts before/after filtering and thresholds used.
    """
    # Work on a copy so the caller's object is never mutated
    adata = adata.copy()

    # Ensure var_names are unique (CITE-seq files mix RNA + ADT)
    adata.var_names_make_unique()

    n_cells_input = adata.n_obs
    n_genes_input = adata.n_vars
    logger.info("QC start — %d cells × %d genes", n_cells_input, n_genes_input)

    # ------------------------------------------------------------------
    # 1. Detect mitochondrial genes
    # ------------------------------------------------------------------
    mt_mask = _detect_mt_genes(adata)
    n_mt = int(mt_mask.sum())
    logger.info("Mitochondrial genes detected: %d", n_mt)
    adata.var["mt"] = mt_mask

    # ------------------------------------------------------------------
    # 2. Compute per-cell QC metrics
    # ------------------------------------------------------------------
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # After this call adata.obs contains:
    #   n_genes_by_counts, total_counts, pct_counts_mt  (and more)

    # ------------------------------------------------------------------
    # 3. Doublet detection with Scrublet
    # ------------------------------------------------------------------
    _run_scrublet(
        adata,
        expected_doublet_rate=scrublet_expected_doublet_rate,
        random_state=random_state,
    )

    # ------------------------------------------------------------------
    # 4. Apply filters
    # ------------------------------------------------------------------
    keep = np.ones(adata.n_obs, dtype=bool)

    low_genes  = np.asarray(adata.obs["n_genes_by_counts"] < min_genes, dtype=bool)
    high_genes = np.asarray(adata.obs["n_genes_by_counts"] > max_genes, dtype=bool)
    high_mt    = np.asarray(adata.obs["pct_counts_mt"] > max_mt_pct, dtype=bool)

    keep &= ~low_genes
    keep &= ~high_genes
    keep &= ~high_mt

    if remove_doublets:
        doublets = np.asarray(adata.obs["predicted_doublet"], dtype=bool)
        keep &= ~doublets
    else:
        doublets = np.zeros(adata.n_obs, dtype=bool)

    adata_filtered = adata[keep].copy()

    n_cells_output = adata_filtered.n_obs
    n_removed = n_cells_input - n_cells_output

    logger.info(
        "QC complete — kept %d / %d cells (removed %d)",
        n_cells_output, n_cells_input, n_removed,
    )

    # ------------------------------------------------------------------
    # 5. Build metrics summary
    # ------------------------------------------------------------------
    metrics: dict[str, Any] = {
        # Counts
        "n_cells_input":       n_cells_input,
        "n_cells_output":      n_cells_output,
        "n_cells_removed":     n_removed,
        "n_genes_input":       n_genes_input,
        # Per-filter removal counts
        "n_removed_low_genes":  int(low_genes.sum()),
        "n_removed_high_genes": int(high_genes.sum()),
        "n_removed_high_mt":    int(high_mt.sum()),
        "n_removed_doublets":   int(doublets.sum()) if remove_doublets else 0,
        # Thresholds used
        "thresholds": {
            "min_genes":       min_genes,
            "max_genes":       max_genes,
            "max_mt_pct":      max_mt_pct,
            "remove_doublets": remove_doublets,
        },
        # Distribution summaries (on pre-filter data)
        "median_genes_per_cell":    float(np.median(adata.obs["n_genes_by_counts"])),
        "median_umi_per_cell":      float(np.median(adata.obs["total_counts"])),
        "median_mt_pct":            float(np.median(adata.obs["pct_counts_mt"])),
        "n_mt_genes":               n_mt,
    }

    # ------------------------------------------------------------------
    # 6. Optional HTML report
    # ------------------------------------------------------------------
    if generate_report:
        try:
            from pipeline.modules.qc.qc_report import generate_qc_report
            generate_qc_report(
                adata_raw=adata,
                adata_filtered=adata_filtered,
                metrics=metrics,
                output_path=report_path,
                sample_name=sample_name,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("QC report generation failed: %s", exc)

    return adata_filtered, metrics


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _detect_mt_genes(adata: AnnData) -> np.ndarray:
    """Return a boolean array (length n_vars) marking mitochondrial genes.

    Checks for the prefixes ``MT-`` (human) and ``mt-`` (mouse).
    Works on both ``adata.var_names`` (gene symbols) and, if present,
    ``adata.var['gene_ids']``.
    """
    names = adata.var_names.astype(str)
    mt_mask = names.str.startswith("MT-") | names.str.startswith("mt-")

    # Fallback: check gene_ids column if no MT genes found in var_names
    if mt_mask.sum() == 0 and "gene_ids" in adata.var.columns:
        gene_ids = adata.var["gene_ids"].astype(str)
        mt_mask = gene_ids.str.startswith("MT-") | gene_ids.str.startswith("mt-")

    # .values is safe on pandas Series; numpy arrays have no .values but are already arrays
    return mt_mask.values if hasattr(mt_mask, 'values') else mt_mask


def _run_scrublet(
    adata: AnnData,
    expected_doublet_rate: float,
    random_state: int,
) -> None:
    """Run Scrublet and attach results to ``adata.obs`` in-place.

    Adds two columns:
      - ``doublet_score``     : continuous score (0–1)
      - ``predicted_doublet`` : boolean flag
    """
    try:
        # Scrublet expects a dense or sparse counts matrix
        counts_matrix = adata.X

        scrub = scr.Scrublet(
            counts_matrix,
            expected_doublet_rate=expected_doublet_rate,
            random_state=random_state,
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

        adata.obs["doublet_score"]     = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets

        n_doublets = int(predicted_doublets.sum())
        logger.info(
            "Scrublet: %d / %d cells flagged as doublets (%.1f%%)",
            n_doublets, adata.n_obs, 100 * n_doublets / adata.n_obs,
        )

    except Exception as exc:  # noqa: BLE001
        # Never let Scrublet failure crash the whole QC run
        logger.warning("Scrublet failed — skipping doublet detection. Reason: %s", exc)
        adata.obs["doublet_score"]     = np.nan
        adata.obs["predicted_doublet"] = False
