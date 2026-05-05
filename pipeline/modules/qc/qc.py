"""
OmicSage — QC Module
pipeline/modules/qc/qc.py

Input:  AnnData with raw counts in adata.X (output of ingest.py).
        May be a plain RNA object, a CITE-seq object (RNA + ADT),
        or a Multiome object (RNA + ATAC).

Output: MuData object + QC metrics dict.

        The MuData always contains a "rna" key.
        CITE-seq additionally contains "adt".
        Multiome additionally contains "atac".

        QC metrics (n_genes_by_counts, total_counts, pct_counts_mt,
        doublet_score, predicted_doublet) are computed on RNA features
        only and live in mdata["rna"].obs.
        Other modalities start clean, ready for their own QC later.

Filters applied (all configurable, all based on RNA metrics):
  - min_genes       : minimum genes per cell (default 200)
  - max_genes       : maximum genes per cell (default 2500)
  - max_mt_pct      : maximum MT% per cell   (default 5.0)
  - remove_doublets : drop Scrublet doublets  (default True)
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import scanpy as sc
import scrublet as scr
from anndata import AnnData
import mudata as md
from mudata import MuData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_qc(
    adata: AnnData,
    modality: str = "auto",
    min_genes: int = 200,
    max_genes: int = 2500,
    max_mt_pct: float = 5.0,
    remove_doublets: bool = True,
    scrublet_expected_doublet_rate: float = 0.06,
    random_state: int = 0,
    generate_report: bool = False,
    report_path: str = "reports/qc_report.html",
    sample_name: str = "sample",
) -> tuple[MuData, dict[str, Any]]:
    """Run QC on a raw-count AnnData object.

    Automatically detects whether the input is plain RNA, CITE-seq
    (RNA + ADT), or Multiome (RNA + ATAC).  QC metrics and cell filters
    are always computed on RNA features only.  The returned MuData
    preserves all modalities with low-quality cells removed.

    Parameters
    ----------
    adata:
        AnnData with raw integer counts in ``adata.X``.
        Produced by ``ingest.py``.
    modality:
        One of ``"auto"`` | ``"rna"`` | ``"cite"`` | ``"multiome"``.
        ``"auto"`` (default) calls ``_detect_modality()`` internally.
    min_genes:
        Minimum number of genes detected per cell (RNA metric).
    max_genes:
        Maximum number of genes detected per cell (RNA metric).
    max_mt_pct:
        Maximum mitochondrial read percentage per cell (RNA metric).
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
    mdata:
        MuData containing one AnnData per modality, keyed by:
        ``"rna"`` (always), ``"adt"`` (CITE-seq), ``"atac"`` (Multiome).
        Only cells that passed all RNA-based QC filters are retained
        across all modalities.
        QC obs columns live on ``mdata["rna"]``.
    metrics:
        Summary dict with counts before/after filtering and thresholds used.
    """
    # Work on a copy so the caller's object is never mutated
    adata = adata.copy()

    # ------------------------------------------------------------------
    # 0. Modality detection
    # ------------------------------------------------------------------
    if modality == "auto":
        modality = _detect_modality(adata)

    logger.info("Detected modality: %s", modality)

    # ------------------------------------------------------------------
    # 1. Split into per-modality AnnData objects
    #    RNA subset is used for all metric computation and filtering.
    #    Other modalities are carried along for the final MuData.
    # ------------------------------------------------------------------
    adata_rna, other_modalities = _split_modalities(adata, modality)

    # Ensure var_names are unique within the RNA subset
    adata_rna.var_names_make_unique()

    n_cells_input = adata_rna.n_obs
    n_genes_input = adata_rna.n_vars
    logger.info("QC start — %d cells × %d RNA features", n_cells_input, n_genes_input)

    # ------------------------------------------------------------------
    # 2. Detect mitochondrial genes (RNA only)
    # ------------------------------------------------------------------
    mt_mask = _detect_mt_genes(adata_rna)
    n_mt = int(mt_mask.sum())
    logger.info("Mitochondrial genes detected: %d", n_mt)
    adata_rna.var["mt"] = mt_mask

    # ------------------------------------------------------------------
    # 3. Compute per-cell QC metrics (RNA only)
    # ------------------------------------------------------------------
    sc.pp.calculate_qc_metrics(
        adata_rna,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # ------------------------------------------------------------------
    # 4. Doublet detection with Scrublet (RNA only)
    # ------------------------------------------------------------------
    _run_scrublet(
        adata_rna,
        expected_doublet_rate=scrublet_expected_doublet_rate,
        random_state=random_state,
    )

    # ------------------------------------------------------------------
    # 5. Build cell keep-mask from RNA metrics
    # ------------------------------------------------------------------
    keep = np.ones(adata_rna.n_obs, dtype=bool)

    low_genes  = np.asarray(adata_rna.obs["n_genes_by_counts"] < min_genes, dtype=bool)
    high_genes = np.asarray(adata_rna.obs["n_genes_by_counts"] > max_genes, dtype=bool)
    high_mt    = np.asarray(adata_rna.obs["pct_counts_mt"] > max_mt_pct, dtype=bool)

    keep &= ~low_genes
    keep &= ~high_genes
    keep &= ~high_mt

    if remove_doublets:
        doublets = np.asarray(adata_rna.obs["predicted_doublet"], dtype=bool)
        keep &= ~doublets
    else:
        doublets = np.zeros(adata_rna.n_obs, dtype=bool)

    keep_barcodes = adata_rna.obs_names[keep]

    # ------------------------------------------------------------------
    # 6. Apply cell filter to RNA and all other modalities
    # ------------------------------------------------------------------
    adata_rna_filtered = adata_rna[keep_barcodes].copy()

    filtered_modalities: dict[str, AnnData] = {"rna": adata_rna_filtered}
    for mod_key, mod_adata in other_modalities.items():
        # Align to the same filtered barcodes
        shared = mod_adata.obs_names.intersection(keep_barcodes)
        filtered_modalities[mod_key] = mod_adata[shared].copy()

    # ------------------------------------------------------------------
    # 7. Assemble MuData
    # ------------------------------------------------------------------
    mdata = MuData(filtered_modalities)

    n_cells_output = adata_rna_filtered.n_obs
    n_removed = n_cells_input - n_cells_output

    logger.info(
        "QC complete — kept %d / %d cells (removed %d)",
        n_cells_output, n_cells_input, n_removed,
    )

    # ------------------------------------------------------------------
    # 8. Build metrics summary
    # ------------------------------------------------------------------
    metrics: dict[str, Any] = {
        # Counts
        "n_cells_input":        n_cells_input,
        "n_cells_output":       n_cells_output,
        "n_cells_removed":      n_removed,
        "n_genes_input":        n_genes_input,
        "modality":             modality,
        # Per-filter removal counts
        "n_removed_low_genes":  int(low_genes.sum()),
        "n_removed_high_genes": int(high_genes.sum()),
        "n_removed_high_mt":    int(high_mt.sum()),
        "n_removed_doublets":   int(doublets.sum()) if remove_doublets else 0,
        # Thresholds used
        "thresholds": {
            "min_genes":        min_genes,
            "max_genes":        max_genes,
            "max_mt_pct":       max_mt_pct,
            "remove_doublets":  remove_doublets,
        },
        # Distribution summaries (on pre-filter RNA data)
        "median_genes_per_cell": float(np.median(adata_rna.obs["n_genes_by_counts"])),
        "median_umi_per_cell":   float(np.median(adata_rna.obs["total_counts"])),
        "median_mt_pct":         float(np.median(adata_rna.obs["pct_counts_mt"])),
        "n_mt_genes":            n_mt,
    }

    # ------------------------------------------------------------------
    # 9. Optional HTML report
    # ------------------------------------------------------------------
    if generate_report:
        try:
            from pipeline.modules.qc.qc_report import generate_qc_report
            generate_qc_report(
                adata_raw=adata_rna,
                adata_filtered=adata_rna_filtered,
                metrics=metrics,
                output_path=report_path,
                sample_name=sample_name,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("QC report generation failed: %s", exc)

    return mdata, metrics


# ---------------------------------------------------------------------------
# Modality detection
# ---------------------------------------------------------------------------

def _detect_modality(adata: AnnData) -> str:
    """Infer whether *adata* is plain RNA, CITE-seq, or Multiome.

    Inspects ``adata.var['feature_types']`` if present.

    Returns
    -------
    ``"rna"``      — no ``feature_types`` column, or all features are RNA.
    ``"cite"``     — ``feature_types`` contains ``"ADT"`` or ``"Antibody Capture"``.
    ``"multiome"`` — ``feature_types`` contains ``"ATAC"`` or ``"Peaks"``.
    """
    if "feature_types" not in adata.var.columns:
        return "rna"

    ft = set(adata.var["feature_types"].astype(str).unique())

    if "ADT" in ft or "Antibody Capture" in ft:
        return "cite"

    if "ATAC" in ft or "Peaks" in ft:
        return "multiome"

    return "rna"


# ---------------------------------------------------------------------------
# Modality splitting
# ---------------------------------------------------------------------------

_GEX_LABELS = {"GEX", "Gene Expression"}
_ADT_LABELS = {"ADT", "Antibody Capture"}
_ATAC_LABELS = {"ATAC", "Peaks"}


def _split_modalities(
    adata: AnnData,
    modality: str,
) -> tuple[AnnData, dict[str, AnnData]]:
    """Split a mixed AnnData into one AnnData per modality.

    Parameters
    ----------
    adata:
        Raw AnnData (may contain multiple feature types).
    modality:
        One of ``"rna"``, ``"cite"``, ``"multiome"``.

    Returns
    -------
    adata_rna:
        AnnData containing only RNA (GEX) features.
    other_modalities:
        Dict mapping modality key → AnnData for non-RNA features.
        Empty for plain RNA data.
    """
    if modality == "rna":
        # Nothing to split — return as-is
        return adata.copy(), {}

    if "feature_types" not in adata.var.columns:
        logger.warning(
            "modality='%s' requested but 'feature_types' column not found in adata.var. "
            "Treating all features as RNA.",
            modality,
        )
        return adata.copy(), {}

    ft = adata.var["feature_types"].astype(str)

    # RNA mask — anything labelled GEX / Gene Expression, or not in a known
    # non-RNA label set (fallback for files with unlabelled RNA rows)
    non_rna_labels = _ADT_LABELS | _ATAC_LABELS
    gex_mask = ft.isin(_GEX_LABELS) | ~ft.isin(non_rna_labels)

    adata_rna = adata[:, gex_mask].copy()
    other: dict[str, AnnData] = {}

    if modality == "cite":
        adt_mask = ft.isin(_ADT_LABELS)
        if adt_mask.sum() == 0:
            logger.warning("CITE-seq modality detected but no ADT features found.")
        else:
            other["adt"] = adata[:, adt_mask].copy()

    elif modality == "multiome":
        atac_mask = ft.isin(_ATAC_LABELS)
        if atac_mask.sum() == 0:
            logger.warning("Multiome modality detected but no ATAC features found.")
        else:
            other["atac"] = adata[:, atac_mask].copy()

    logger.info(
        "Split: %d RNA features, %s",
        adata_rna.n_vars,
        ", ".join(f"{v.n_vars} {k.upper()} features" for k, v in other.items())
        or "no additional modalities",
    )

    return adata_rna, other


# ---------------------------------------------------------------------------
# Private helpers (unchanged from v1)
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

    return mt_mask.values if hasattr(mt_mask, "values") else mt_mask


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
        scrub = scr.Scrublet(
            adata.X,
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
        logger.warning("Scrublet failed — skipping doublet detection. Reason: %s", exc)
        adata.obs["doublet_score"]     = np.nan
        adata.obs["predicted_doublet"] = False
