"""
OmicSage — ATAC QC Module
pipeline/modules/multiome/atac_qc.py

Input : AnnData with raw ATAC peak counts in .X  (mdata["atac"] from run_qc())
        For GSE194122 multiome this is the NeurIPS 2021 BMMC dataset where
        CellRanger-ARC pre-computed fragment-derived metrics already live in .obs.

Output: AnnData + QC metrics dict

What this module does
---------------------
1. Rename CellRanger-ARC obs columns to OmicSage-namespaced equivalents
   (ATAC_nCount_peaks → total_peak_counts, etc.)
2. Preserve ground-truth labels as obs["cell_type_groundtruth"]
3. Compute n_peaks_by_counts via sc.pp.calculate_qc_metrics (not in CellRanger obs)
4. Run Scrublet on the peak matrix → atac_predicted_doublet (flag only)
5. Apply QC filters (flag only by default — RNA QC already filtered cells)
6. Filter features: keep peaks present in ≥ min_cells cells
7. Save raw counts to .layers["counts"]
8. Write provenance to .uns["omicsage_atac_qc"]

CellRanger-ARC obs columns used (when present)
-----------------------------------------------
    ATAC_nCount_peaks        → total_peak_counts
    ATAC_atac_fragments      → total_fragment_counts
    ATAC_reads_in_peaks_frac → reads_in_peaks_frac
    ATAC_blacklist_fraction  → blacklist_fraction
    ATAC_nucleosome_signal   → nucleosome_signal

TSS enrichment and nucleosome signal from raw fragment files are NOT computed
here — the GSE194122 multiome data is GEO-deposited processed data without
fragment files.  CellRanger-ARC nucleosome_signal values are re-used directly.

Thresholds
----------
Default thresholds are calibrated for the NeurIPS 2021 BMMC multiome dataset,
taken directly from the sc-best-practices ATAC QC chapter (Lance & Martens 2022)
which uses this exact dataset:
    total_peak_counts :   > 1500  and < 100 000
    n_peaks_by_counts :   > 750
    nucleosome_signal :   < 2
    (TSS enrichment not computed — no fragment file)

Usage
-----
    from pipeline.modules.multiome.atac_qc import atac_qc

    atac = sc.read_h5ad("data/processed/GSE194122_multiome/01_qc_atac.h5ad")
    atac_qcd, metrics = atac_qc(atac)
"""

from __future__ import annotations

import logging
import warnings
from datetime import datetime, timezone
from typing import Any, Optional

import numpy as np
import scanpy as sc
import scrublet as scr
import scipy.sparse as sp
from anndata import AnnData

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CellRanger-ARC column → OmicSage column mapping
# ---------------------------------------------------------------------------
_CELLRANGER_COL_MAP = {
    "ATAC_nCount_peaks":        "total_peak_counts",
    "ATAC_atac_fragments":      "total_fragment_counts",
    "ATAC_reads_in_peaks_frac": "reads_in_peaks_frac",
    "ATAC_blacklist_fraction":  "blacklist_fraction",
    "ATAC_nucleosome_signal":   "nucleosome_signal",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def atac_qc(
    adata: AnnData,
    # Cell-level thresholds (sc-best-practices NeurIPS 2021 BMMC values)
    min_peaks: int = 750,
    max_peaks: int = 500_000,
    min_peak_counts: int = 1_500,
    max_peak_counts: int = 100_000,
    max_nucleosome_signal: float = 2.0,
    # Feature-level threshold
    min_cells: int = 15,
    # Doublet detection
    run_scrublet: bool = True,
    scrublet_expected_doublet_rate: float = 0.06,
    random_state: int = 0,
    # Behaviour
    filter_cells: bool = False,
    inplace: bool = False,
) -> tuple[AnnData, dict[str, Any]]:
    """Run ATAC-specific QC on a peak count AnnData.

    Designed for the output of ``run_qc()`` (``mdata["atac"]``) for multiome
    datasets.  Cell filtering is disabled by default because the RNA QC step
    has already removed low-quality barcodes; this step computes and records
    ATAC-specific metrics only.

    Parameters
    ----------
    adata : AnnData
        Raw ATAC AnnData with peak counts in ``.X``.
        Typically ``mdata["atac"]`` loaded from ``01_qc_atac.h5ad``.
        CellRanger-ARC pre-computed metrics are expected in ``.obs`` (see
        ``_CELLRANGER_COL_MAP``).
    min_peaks : int
        Minimum number of peaks with > 0 counts per cell.
    max_peaks : int
        Maximum number of peaks with > 0 counts per cell.
    min_peak_counts : int
        Minimum total peak counts per cell (depth filter).
    max_peak_counts : int
        Maximum total peak counts per cell (doublet / artifact filter).
    max_nucleosome_signal : float
        Maximum nucleosome signal per cell.  Values > threshold indicate
        poor signal-to-noise ratio.  Only applied when
        ``nucleosome_signal`` is available (CellRanger-ARC or fragment file).
    min_cells : int
        Minimum number of cells a peak must be detected in to be retained.
        sc-best-practices threshold: 15 cells (derived from 2% rare cell type
        at 5% detection rate in 15 000 cells).
    run_scrublet : bool
        If True, run Scrublet on the peak count matrix and add
        ``atac_predicted_doublet`` and ``atac_doublet_score`` to ``.obs``.
        Flag only — cells are never removed based on this score here.
    scrublet_expected_doublet_rate : float
        Prior doublet rate passed to Scrublet.
    random_state : int
        Reproducibility seed.
    filter_cells : bool
        If True, remove cells that fail QC thresholds.
        Default False — RNA QC already filtered cells; this step flags only.
    inplace : bool
        If False (default), operate on a copy.

    Returns
    -------
    adata_out : AnnData
        QC-annotated ATAC AnnData.
        New ``.obs`` columns (always):
          - ``n_peaks_by_counts``     number of peaks > 0 per cell
          - ``total_peak_counts``     total counts per cell
          - ``cell_type_groundtruth`` preserved from ``cell_type`` if present
        CellRanger-ARC columns renamed to OmicSage names (when present):
          - ``total_fragment_counts``, ``reads_in_peaks_frac``,
            ``blacklist_fraction``, ``nucleosome_signal``
        QC flag columns:
          - ``atac_qc_pass``           True if cell passes all filters
          - ``atac_doublet_score``      Scrublet score (if run_scrublet=True)
          - ``atac_predicted_doublet``  Scrublet flag  (if run_scrublet=True)
        ``.layers["counts"]``  raw peak counts (always written)
        ``.uns["omicsage_atac_qc"]``  provenance record
    metrics : dict
        Summary dict with cell/peak counts before/after and thresholds used.

    Raises
    ------
    TypeError
        If ``adata`` is not an AnnData object.
    ValueError
        If ``adata`` has 0 cells or 0 peaks.
    """
    _validate_input(adata)

    # Always work on a concrete (non-view) object so that uns writes are
    # not silently dropped by anndata's copy-on-write mechanism.
    adata_out = adata.copy()

    n_cells_before = adata_out.n_obs
    n_peaks_before = adata_out.n_vars
    logger.info(
        "ATAC QC start — %d cells × %d peaks", n_cells_before, n_peaks_before
    )

    # ------------------------------------------------------------------
    # Step 1 — Preserve ground-truth labels
    # ------------------------------------------------------------------
    if "cell_type" in adata_out.obs.columns:
        adata_out.obs["cell_type_groundtruth"] = adata_out.obs["cell_type"].copy()
        logger.info("Preserved obs['cell_type'] → obs['cell_type_groundtruth']")

    # ------------------------------------------------------------------
    # Step 2 — Rename CellRanger-ARC columns to OmicSage names
    # ------------------------------------------------------------------
    cellranger_cols_found: list[str] = []
    for cr_col, omicsage_col in _CELLRANGER_COL_MAP.items():
        if cr_col in adata_out.obs.columns:
            adata_out.obs[omicsage_col] = adata_out.obs[cr_col].copy()
            cellranger_cols_found.append(cr_col)
            logger.debug("Renamed obs['%s'] → obs['%s']", cr_col, omicsage_col)

    cellranger_metrics_available = len(cellranger_cols_found) > 0
    logger.info(
        "CellRanger-ARC metrics found: %s",
        cellranger_cols_found if cellranger_metrics_available else "none",
    )

    # ------------------------------------------------------------------
    # Step 3 — Compute n_peaks_by_counts via scanpy
    #          (not the same as ATAC_nCount_peaks — n_peaks counts
    #           peaks with > 0 entries, nCount_peaks sums all counts)
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.calculate_qc_metrics(
            adata_out, percent_top=None, log1p=False, inplace=True
        )

    # scanpy names these for RNA — rename to ATAC convention
    if "n_genes_by_counts" in adata_out.obs.columns:
        adata_out.obs["n_peaks_by_counts"] = adata_out.obs["n_genes_by_counts"].copy()
        del adata_out.obs["n_genes_by_counts"]

    # If CellRanger didn't provide total_peak_counts, use scanpy's total_counts
    if "total_peak_counts" not in adata_out.obs.columns:
        if "total_counts" in adata_out.obs.columns:
            adata_out.obs["total_peak_counts"] = adata_out.obs["total_counts"].copy()

    logger.info(
        "Computed n_peaks_by_counts — median=%.0f",
        float(np.median(adata_out.obs["n_peaks_by_counts"])),
    )

    # ------------------------------------------------------------------
    # Step 4 — Scrublet doublet detection on peak matrix
    # ------------------------------------------------------------------
    if run_scrublet:
        _run_scrublet(
            adata_out,
            expected_doublet_rate=scrublet_expected_doublet_rate,
            random_state=random_state,
        )
    else:
        adata_out.obs["atac_doublet_score"]     = np.nan
        adata_out.obs["atac_predicted_doublet"] = False

    # ------------------------------------------------------------------
    # Step 5 — Build per-cell QC pass/fail mask
    # ------------------------------------------------------------------
    keep = np.ones(adata_out.n_obs, dtype=bool)

    low_peaks  = np.asarray(adata_out.obs["n_peaks_by_counts"] < min_peaks,       dtype=bool)
    high_peaks = np.asarray(adata_out.obs["n_peaks_by_counts"] > max_peaks,       dtype=bool)
    low_counts = np.asarray(adata_out.obs["total_peak_counts"] < min_peak_counts, dtype=bool)
    high_counts= np.asarray(adata_out.obs["total_peak_counts"] > max_peak_counts, dtype=bool)

    keep &= ~low_peaks
    keep &= ~high_peaks
    keep &= ~low_counts
    keep &= ~high_counts

    # Nucleosome signal filter (only when column is available)
    high_nuc = np.zeros(adata_out.n_obs, dtype=bool)
    nucleosome_filter_applied = False
    if "nucleosome_signal" in adata_out.obs.columns:
        high_nuc = np.asarray(
            adata_out.obs["nucleosome_signal"] > max_nucleosome_signal, dtype=bool
        )
        keep &= ~high_nuc
        nucleosome_filter_applied = True

    adata_out.obs["atac_qc_pass"] = keep

    n_pass = int(keep.sum())
    n_fail = adata_out.n_obs - n_pass
    logger.info(
        "QC flags — %d pass, %d fail (filter_cells=%s)",
        n_pass, n_fail, filter_cells,
    )

    # ------------------------------------------------------------------
    # Step 6 — Optionally remove failing cells
    # ------------------------------------------------------------------
    if filter_cells and n_fail > 0:
        adata_out = adata_out[keep].copy()
        logger.info("Removed %d low-quality ATAC cells", n_fail)

    n_cells_after_cell_filter = adata_out.n_obs

    # ------------------------------------------------------------------
    # Step 7 — Feature filter: keep peaks in ≥ min_cells cells
    # ------------------------------------------------------------------
    if "n_cells_by_counts" in adata_out.var.columns:
        peak_keep = adata_out.var["n_cells_by_counts"] >= min_cells
    else:
        # fallback: compute from X
        if sp.issparse(adata_out.X):
            cells_per_peak = np.asarray((adata_out.X > 0).sum(axis=0)).ravel()
        else:
            cells_per_peak = (adata_out.X > 0).sum(axis=0)
        peak_keep = cells_per_peak >= min_cells

    n_peaks_removed = int((~peak_keep).sum())
    adata_out = adata_out[:, peak_keep].copy()
    n_peaks_after = adata_out.n_vars
    logger.info(
        "Feature filter — removed %d peaks present in < %d cells  (%d remain)",
        n_peaks_removed, min_cells, n_peaks_after,
    )

    # ------------------------------------------------------------------
    # Step 8 — Save raw counts as a layer
    # ------------------------------------------------------------------
    if "counts" not in adata_out.layers:
        adata_out.layers["counts"] = adata_out.X.copy()
        logger.info("Raw peak counts saved to adata.layers['counts']")

    # ------------------------------------------------------------------
    # Step 9 — Build metrics dict
    # ------------------------------------------------------------------
    metrics: dict[str, Any] = {
        "n_cells_before":               n_cells_before,
        "n_cells_after":                n_cells_after_cell_filter,
        "n_cells_removed_qc":           n_cells_before - n_cells_after_cell_filter,
        "n_peaks_before":               n_peaks_before,
        "n_peaks_after":                n_peaks_after,
        "n_peaks_removed_feature_filter": n_peaks_removed,
        "n_qc_pass":                    n_pass,
        "n_qc_fail":                    n_fail,
        "filter_cells_applied":         filter_cells,
        "cellranger_metrics_available": cellranger_metrics_available,
        "cellranger_cols_found":        cellranger_cols_found,
        "fragment_file_available":      False,
        "nucleosome_filter_applied":    nucleosome_filter_applied,
        # Per-filter breakdown
        "n_fail_low_peaks":             int(low_peaks.sum()),
        "n_fail_high_peaks":            int(high_peaks.sum()),
        "n_fail_low_counts":            int(low_counts.sum()),
        "n_fail_high_counts":           int(high_counts.sum()),
        "n_fail_high_nucleosome":       int(high_nuc.sum()),
        # Distribution summaries (post-rename, pre-filter)
        "median_n_peaks_by_counts":     float(np.median(
            adata_out.obs["n_peaks_by_counts"])) if n_cells_after_cell_filter > 0 else 0.0,
        "median_total_peak_counts":     float(np.median(
            adata_out.obs["total_peak_counts"])) if n_cells_after_cell_filter > 0 else 0.0,
        # Thresholds used
        "thresholds": {
            "min_peaks":              min_peaks,
            "max_peaks":              max_peaks,
            "min_peak_counts":        min_peak_counts,
            "max_peak_counts":        max_peak_counts,
            "max_nucleosome_signal":  max_nucleosome_signal,
            "min_cells_per_peak":     min_cells,
        },
    }

    # nucleosome_signal distribution summary (when available)
    if "nucleosome_signal" in adata_out.obs.columns:
        metrics["median_nucleosome_signal"] = float(
            np.median(adata_out.obs["nucleosome_signal"])
        )

    # ------------------------------------------------------------------
    # Step 10 — Provenance
    # ------------------------------------------------------------------
    adata_out.uns["omicsage_atac_qc"] = {
        "omicsage_module":              "pipeline.modules.multiome.atac_qc",
        "omicsage_version":             "0.1.0",
        "timestamp":                    datetime.now(timezone.utc).isoformat(),
        "n_cells_before":               n_cells_before,
        "n_cells_after":                n_cells_after_cell_filter,
        "n_peaks_before":               n_peaks_before,
        "n_peaks_after":                n_peaks_after,
        "fragment_file_available":      False,
        "cellranger_metrics_used":      cellranger_metrics_available,
        "cellranger_cols_renamed":      _CELLRANGER_COL_MAP,
        "nucleosome_filter_applied":    nucleosome_filter_applied,
        "filter_cells":                 filter_cells,
        "thresholds":                   metrics["thresholds"],
        "note": (
            "TSS enrichment and nucleosome signal from fragment files not computed. "
            "CellRanger-ARC pre-computed metrics used where available. "
            "Cell filtering skipped — RNA QC already removed low-quality barcodes."
            if not filter_cells else
            "Cell filtering applied based on ATAC-specific QC metrics."
        ),
    }

    logger.info(
        "ATAC QC complete — %d cells, %d peaks, %d QC pass (%.1f%%)",
        n_cells_after_cell_filter,
        n_peaks_after,
        n_pass,
        100 * n_pass / n_cells_before if n_cells_before > 0 else 0.0,
    )

    # When inplace=True, write all results back to the caller's object.
    # We always work on a copy internally (to avoid view/uns issues), then
    # copy results back so the caller's reference reflects all changes.
    if inplace:
        for col in adata_out.obs.columns:
            adata.obs[col] = adata_out.obs[col]
        for key in adata_out.uns:
            adata.uns[key] = adata_out.uns[key]
        for key in adata_out.layers:
            if adata_out.n_vars == adata.n_vars:   # only if peak set unchanged
                adata.layers[key] = adata_out.layers[key]

    return adata_out, metrics


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _validate_input(adata: AnnData) -> None:
    if not isinstance(adata, AnnData):
        raise TypeError(
            f"atac_qc() expects an AnnData object, got {type(adata).__name__}. "
            "Pass mdata['atac'] not the full MuData."
        )
    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 cells — nothing to QC.")
    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 peaks — nothing to QC.")


def _run_scrublet(
    adata: AnnData,
    expected_doublet_rate: float,
    random_state: int,
) -> None:
    """Run Scrublet on peak matrix and write results to adata.obs in-place."""
    try:
        # Scrublet requires dense or CSR; ensure CSR
        X = adata.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        elif not sp.isspmatrix_csr(X):
            X = X.tocsr()

        scrub = scr.Scrublet(
            X,
            expected_doublet_rate=expected_doublet_rate,
            random_state=random_state,
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

        adata.obs["atac_doublet_score"]     = doublet_scores
        adata.obs["atac_predicted_doublet"] = predicted_doublets

        n_doublets = int(predicted_doublets.sum())
        logger.info(
            "Scrublet (ATAC): %d / %d cells flagged as doublets (%.1f%%)",
            n_doublets, adata.n_obs, 100 * n_doublets / adata.n_obs,
        )

    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "Scrublet failed on ATAC peak matrix — skipping. Reason: %s", exc
        )
        adata.obs["atac_doublet_score"]     = np.nan
        adata.obs["atac_predicted_doublet"] = False
