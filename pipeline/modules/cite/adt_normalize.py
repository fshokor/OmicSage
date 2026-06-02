"""
OmicSage — ADT Normalization Module
pipeline/modules/cite/adt_normalize.py

Normalizes raw ADT (antibody-derived tag) count data from CITE-seq experiments.
Input : AnnData with raw integer ADT counts in .X  (typically mdata["adt"])
Output: updated AnnData  +  metrics dict

Methods supported
-----------------
CLR  (Centered Log-Ratio) — always applied
     Divides each count by the geometric mean across all proteins in the same
     axis, then applies log1p.  Two axis conventions are supported:
       axis=0  CLR across cells per protein  (muon default, Seurat margin=1)
       axis=1  CLR across proteins per cell  (Seurat margin=2)
     The sc-best-practices reference calls mu.prot.pp.clr without an axis
     argument, so axis=0 is the default here too.

DSB  (Denoised and Scaled by Background) — optional, requires empty droplets
     Removes ambient noise using empty droplets and unspecific binding noise
     using isotype controls.  Pass dsb_empty_adata to enable.
     Implemented via muon.prot.pp.dsb (Mulè et al. 2022, Nat Comm 13:2099).
     After DSB, adata.X is set to the DSB-normalized values so all downstream
     steps (PCA, doublet detection, annotation) operate on DSB by default.

Usage
-----
    from pipeline.modules.cite.adt_normalize import normalize_adt

    # CLR only (no empty droplets available):
    adata_norm, metrics = normalize_adt(mdata["adt"])

    # DSB + CLR (recommended when empty droplets are available):
    adata_norm, metrics = normalize_adt(
        mdata["adt"],
        dsb_empty_adata=empty_adt,          # AnnData of empty droplets
        isotype_controls=["IgG1", "IgG2a"], # optional
    )

Layers written
--------------
    .layers["counts"]   — original raw integer counts (always)
    .layers["adt_clr"]  — CLR-normalised values (always)
    .layers["adt_dsb"]  — DSB-normalised values (only when dsb_empty_adata provided)
    .X                  — DSB values if DSB was run, else CLR values

Provenance
----------
    .uns["omicsage_adt_normalize"] — full parameter + metrics record
"""

from __future__ import annotations

import logging
import warnings
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import muon
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def normalize_adt(
    adata: AnnData,
    clr_axis: int = 0,
    dsb_empty_adata: Optional[AnnData] = None,
    isotype_controls: Optional[list[str]] = None,
    dsb_pseudocount: int = 10,
    dsb_denoise: bool = True,
    dsb_random_state: int = 0,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Normalize raw ADT counts from a CITE-seq experiment.

    CLR normalization is always applied.  DSB normalization is applied in
    addition when ``dsb_empty_adata`` is provided, and its output becomes
    the default ``adata.X`` for downstream steps.

    Parameters
    ----------
    adata : AnnData
        AnnData with raw integer ADT counts in ``.X``.
        Typically ``mdata["adt"]`` from ``run_qc()``.
    clr_axis : int
        Axis along which CLR is computed.
        ``0``  (default) — CLR across cells per protein  (muon default).
        ``1``  — CLR across proteins per cell.
    dsb_empty_adata : AnnData, optional
        AnnData containing ADT counts from **empty droplets** (the unfiltered
        droplet pool minus cell-containing barcodes).  Must have the same
        proteins (var_names) as ``adata``.  When ``None`` (default), DSB is
        skipped and ``adata.X`` is set to the CLR layer.
    isotype_controls : list of str, optional
        Names of isotype control antibodies present in ``adata.var_names``.
        Used by DSB to model and remove unspecific binding noise.
        Ignored when ``dsb_empty_adata`` is ``None``.
    dsb_pseudocount : int
        Pseudocount added before log-transform inside DSB (default 10,
        matching the DSB paper recommendation).
    dsb_denoise : bool
        Whether to apply the per-cell denoising step inside DSB (default True).
        Set to False to speed up testing on small fixtures.
    dsb_random_state : int
        Random seed for DSB Gaussian mixture models (default 0).
    inplace : bool
        If ``False`` (default) work on a copy so the caller's object is not
        mutated.

    Returns
    -------
    adata_out : AnnData
        Normalized ADT AnnData.

        - ``.X``                            — DSB values (if DSB run) else CLR
        - ``.layers["counts"]``             — original raw integer counts
        - ``.layers["adt_clr"]``            — CLR-normalised values
        - ``.layers["adt_dsb"]``            — DSB-normalised values (if run)
        - ``.uns["omicsage_adt_normalize"]`` — provenance record
    metrics : dict
        Summary statistics suitable for logging and downstream QA.

    Raises
    ------
    TypeError
        If ``adata`` is not an AnnData object.
    ValueError
        If ``adata`` has 0 cells or 0 features, or if ``clr_axis`` is not
        0 or 1, or if ``dsb_empty_adata`` has different proteins than ``adata``.
    """
    _validate_input(adata, clr_axis, dsb_empty_adata)

    adata_out = adata if inplace else adata.copy()

    # ------------------------------------------------------------------
    # Step 1 — Preserve raw counts before any modification
    # ------------------------------------------------------------------
    if "counts" not in adata_out.layers:
        logger.info("Saving raw ADT counts to adata.layers['counts']")
        adata_out.layers["counts"] = adata_out.X.copy()
    else:
        logger.info("adata.layers['counts'] already exists — skipping copy")

    # ------------------------------------------------------------------
    # Step 2 — CLR normalization (always run, stored as fallback)
    # ------------------------------------------------------------------
    logger.info(
        "Applying CLR normalization (axis=%d, %s)",
        clr_axis,
        "per-protein across cells" if clr_axis == 0 else "per-cell across proteins",
    )

    # CLR modifies .X in-place; restore raw counts first so DSB (if run)
    # also starts from raw integer counts
    raw_copy = adata_out.layers["counts"].copy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        muon.prot.pp.clr(adata_out, axis=clr_axis, inplace=True)

    adata_out.layers["adt_clr"] = adata_out.X.copy()
    logger.info("CLR values saved to adata.layers['adt_clr']")

    # Restore raw counts to .X so DSB operates on integer counts
    adata_out.X = raw_copy

    # ------------------------------------------------------------------
    # Step 3 — DSB normalization (optional)
    # ------------------------------------------------------------------
    dsb_applied = False
    dsb_metrics: dict = {}

    if dsb_empty_adata is not None:
        logger.info(
            "Applying DSB normalization (%d empty droplets, denoise=%s)",
            dsb_empty_adata.n_obs,
            dsb_denoise,
        )

        # Filter isotype controls to only those present in var_names
        iso_controls_used: Optional[list[str]] = None
        if isotype_controls:
            present = [c for c in isotype_controls if c in adata_out.var_names]
            missing = [c for c in isotype_controls if c not in adata_out.var_names]
            if missing:
                logger.warning(
                    "Isotype controls not found in var_names (will be ignored): %s",
                    missing,
                )
            iso_controls_used = present if present else None

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            muon.prot.pp.dsb(
                data=adata_out,
                data_raw=dsb_empty_adata,
                pseudocount=dsb_pseudocount,
                denoise_counts=dsb_denoise,
                isotype_controls=iso_controls_used,
                add_layer=True,       # writes to .layers["dsb"], not .X
                random_state=dsb_random_state,
            )

        # muon writes to .layers["dsb"] — rename to OmicSage convention
        if "dsb" in adata_out.layers:
            adata_out.layers["adt_dsb"] = adata_out.layers.pop("dsb")
        else:
            raise RuntimeError(
                "muon.prot.pp.dsb completed but did not write .layers['dsb']. "
                "Check muon version compatibility."
            )

        # Set .X to DSB values — this is what all downstream steps will use
        adata_out.X = adata_out.layers["adt_dsb"].copy()
        logger.info("DSB values saved to adata.layers['adt_dsb'] and set as adata.X")

        dsb_applied = True
        dsb_arr = (
            adata_out.layers["adt_dsb"].toarray()
            if hasattr(adata_out.layers["adt_dsb"], "toarray")
            else np.asarray(adata_out.layers["adt_dsb"])
        )
        dsb_metrics = {
            "dsb_n_empty_droplets": int(dsb_empty_adata.n_obs),
            "dsb_pseudocount": dsb_pseudocount,
            "dsb_denoise": dsb_denoise,
            "dsb_isotype_controls_used": iso_controls_used or [],
            "dsb_mean": float(dsb_arr.mean()),
            "dsb_max": float(dsb_arr.max()),
            "dsb_min": float(dsb_arr.min()),
            # Fraction of protein-cell combinations that are DSB-positive
            # (> 0 = above background threshold, i.e. "truly expressed")
            "dsb_fraction_positive": float((dsb_arr > 0).mean()),
        }
        logger.info(
            "DSB complete — %.1f%% protein-cell combinations above background",
            dsb_metrics["dsb_fraction_positive"] * 100,
        )
    else:
        # No DSB — set .X to CLR (pipeline default when no empty droplets)
        adata_out.X = adata_out.layers["adt_clr"].copy()
        logger.info("DSB skipped — adata.X set to CLR values")

    # ------------------------------------------------------------------
    # Step 4 — Build metrics dict
    # ------------------------------------------------------------------
    x_arr = (
        adata_out.X.toarray()
        if hasattr(adata_out.X, "toarray")
        else np.asarray(adata_out.X)
    )
    raw_arr = (
        adata_out.layers["counts"].toarray()
        if hasattr(adata_out.layers["counts"], "toarray")
        else np.asarray(adata_out.layers["counts"])
    )
    clr_arr = (
        adata_out.layers["adt_clr"].toarray()
        if hasattr(adata_out.layers["adt_clr"], "toarray")
        else np.asarray(adata_out.layers["adt_clr"])
    )

    metrics: dict = {
        "n_cells": int(adata_out.n_obs),
        "n_proteins": int(adata_out.n_vars),
        "clr_axis": clr_axis,
        "dsb_applied": dsb_applied,
        "active_layer": "adt_dsb" if dsb_applied else "adt_clr",
        "raw_counts_in_layer": "counts",
        "clr_in_layer": "adt_clr",
        # CLR distribution summaries
        "clr_mean_per_cell": float(clr_arr.sum(axis=1).mean()),
        "clr_max": float(clr_arr.max()),
        "clr_min": float(clr_arr.min()),
        # Raw count summaries
        "raw_median_total_counts_per_cell": float(
            np.median(raw_arr.sum(axis=1))
        ),
        "raw_max_count": float(raw_arr.max()),
    }
    if dsb_applied:
        metrics.update(dsb_metrics)

    # ------------------------------------------------------------------
    # Step 5 — Store provenance in uns
    # ------------------------------------------------------------------
    adata_out.uns["omicsage_adt_normalize"] = {
        "clr_axis": clr_axis,
        "clr_axis_description": (
            "per-protein across cells (muon default)"
            if clr_axis == 0
            else "per-cell across proteins"
        ),
        "dsb_applied": dsb_applied,
        "dsb_params": dsb_metrics if dsb_applied else {},
        "active_layer": "adt_dsb" if dsb_applied else "adt_clr",
        "normalized_in_layer": "adt_dsb" if dsb_applied else "adt_clr",
        "raw_in_layer": "counts",
        "muon_version": muon.__version__,
        "omicsage_module": "pipeline.modules.cite.adt_normalize",
        "omicsage_version": "0.1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    logger.info(
        "ADT normalization complete — %d cells, %d proteins, CLR axis=%d, DSB=%s",
        adata_out.n_obs,
        adata_out.n_vars,
        clr_axis,
        dsb_applied,
    )

    return adata_out, metrics


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_input(
    adata: AnnData,
    clr_axis: int,
    dsb_empty_adata: Optional[AnnData],
) -> None:
    """Raise informative errors for common input mistakes."""
    if not isinstance(adata, AnnData):
        raise TypeError(
            f"normalize_adt() expects an AnnData object, got {type(adata).__name__}. "
            "Pass mdata['adt'] not the full MuData."
        )

    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 cells — nothing to normalize.")

    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 proteins — nothing to normalize.")

    if clr_axis not in (0, 1):
        raise ValueError(
            f"clr_axis must be 0 or 1, got {clr_axis!r}. "
            "Use 0 for per-protein CLR (muon default) or 1 for per-cell CLR."
        )

    if dsb_empty_adata is not None:
        if not isinstance(dsb_empty_adata, AnnData):
            raise TypeError(
                f"dsb_empty_adata must be an AnnData object, got "
                f"{type(dsb_empty_adata).__name__}."
            )
        if dsb_empty_adata.n_obs == 0:
            raise ValueError(
                "dsb_empty_adata has 0 rows — need at least a few hundred "
                "empty droplets for DSB to estimate ambient background reliably. "
                "Pass dsb_empty_adata=None to fall back to CLR only."
            )
        if dsb_empty_adata.n_vars != adata.n_vars:
            raise ValueError(
                f"dsb_empty_adata has {dsb_empty_adata.n_vars} proteins but "
                f"adata has {adata.n_vars}. Both must have the same proteins."
            )
        # Check var_names match (order matters for muon DSB)
        if list(dsb_empty_adata.var_names) != list(adata.var_names):
            raise ValueError(
                "dsb_empty_adata.var_names do not match adata.var_names. "
                "Subset and reorder dsb_empty_adata to match before calling normalize_adt()."
            )

    # Warn (don't raise) if .X does not look like raw integer counts
    x_sample = adata.X[:min(100, adata.n_obs), :min(100, adata.n_vars)]
    if hasattr(x_sample, "toarray"):
        x_sample = x_sample.toarray()
    x_sample = np.asarray(x_sample, dtype=float)

    non_integer_frac = np.mean(np.abs(x_sample - np.round(x_sample)) > 1e-3)
    if non_integer_frac > 0.05:
        logger.warning(
            "%.0f%% of sampled ADT values are non-integer. "
            "CLR and DSB normalization expect raw counts in .X. "
            "If data is already normalized, results may be unreliable.",
            non_integer_frac * 100,
        )

    logger.debug(
        "Input validation passed — %d cells × %d proteins", adata.n_obs, adata.n_vars
    )
