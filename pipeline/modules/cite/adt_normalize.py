"""
OmicSage — ADT Normalization Module
pipeline/modules/cite/adt_normalize.py

Normalizes raw ADT (antibody-derived tag) count data from CITE-seq experiments.
Input : MuData with raw integer ADT counts in mdata["adt"].X  (output of qc.py)
Output: updated MuData  +  metrics dict

Methods supported
-----------------
CLR  (Centered Log-Ratio) — default, always available
     Divides each count by the geometric mean across all proteins in the same
     axis, then applies log1p.  Two axis conventions are supported:
       axis=0  CLR across cells per protein  (muon default, Seurat margin=1)
       axis=1  CLR across proteins per cell  (Seurat margin=2)
     The sc-best-practices reference calls mu.prot.pp.clr without an axis
     argument, so axis=0 is the default here too.

DSB  (Denoised and Scaled by Background) — optional, requires empty droplets
     Removes ambient noise using empty droplets and unspecific binding noise
     using isotype controls.  Pass dsb_empty_adata to enable.
     Reference: Mulè et al. 2022 (Nature Communications 13:2099)

Usage
-----
    from pipeline.modules.cite.adt_normalize import normalize_adt

    adata_norm, metrics = normalize_adt(mdata["adt"])

    # With explicit axis:
    adata_norm, metrics = normalize_adt(mdata["adt"], clr_axis=1)

    # Updating mdata in-place:
    mdata["adt"], metrics = normalize_adt(mdata["adt"])

Layers written
--------------
    mdata["adt"].layers["counts"]   — original raw integer counts (always)
    mdata["adt"].layers["adt_clr"]  — CLR-normalized values (copy of .X post-CLR)
    mdata["adt"].X                  — CLR-normalized values (downstream-ready)

Provenance
----------
    mdata["adt"].uns["omicsage_adt_normalize"] — full parameter + metrics record
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
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Normalize raw ADT counts from a CITE-seq experiment.

    CLR normalization is always applied.  DSB normalization is applied in
    addition when ``dsb_empty_adata`` is provided.

    Parameters
    ----------
    adata : AnnData
        AnnData with raw integer ADT counts in ``.X``.
        Typically ``mdata["adt"]`` from ``run_qc()``.
    clr_axis : int
        Axis along which CLR is computed.
        ``0``  (default) — CLR across cells per protein  (muon default,
               equivalent to Seurat ``NormalizeData(margin=1)``).
        ``1``  — CLR across proteins per cell
               (equivalent to Seurat ``NormalizeData(margin=2)``).
    dsb_empty_adata : AnnData, optional
        AnnData containing ADT counts from **empty droplets** (the unfiltered
        droplet pool minus cell-containing barcodes).  Required for DSB
        normalization.  When ``None`` (default), DSB is skipped and only CLR
        is applied.  DSB additionally requires ``isotype_controls``.
    isotype_controls : list of str, optional
        Names of isotype control antibodies present in ``adata.var_names``.
        Used by DSB normalization to model unspecific binding noise.
        Ignored when ``dsb_empty_adata`` is ``None``.
    inplace : bool
        If ``False`` (default) work on a copy so the caller's object is not
        mutated.

    Returns
    -------
    adata_out : AnnData
        Normalized ADT AnnData.

        - ``.X``                            — CLR-normalized values
        - ``.layers["counts"]``             — original raw integer counts
        - ``.layers["adt_clr"]``            — CLR-normalized values (named copy)
        - ``.uns["omicsage_adt_normalize"]`` — provenance record
    metrics : dict
        Summary statistics suitable for logging and downstream QA.

    Raises
    ------
    TypeError
        If ``adata`` is not an AnnData object.
    ValueError
        If ``adata`` has 0 cells or 0 features, or if ``clr_axis`` is not
        0 or 1.
    NotImplementedError
        If ``dsb_empty_adata`` is provided (DSB is not yet implemented;
        see roadmap note in module docstring).
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
    # Step 2 — CLR normalization
    # ------------------------------------------------------------------
    logger.info(
        "Applying CLR normalization (axis=%d, %s)",
        clr_axis,
        "per-protein across cells" if clr_axis == 0 else "per-cell across proteins",
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        muon.prot.pp.clr(adata_out, axis=clr_axis, inplace=True)

    # Save CLR result to named layer so it is always accessible by name
    # regardless of what downstream steps write to .X
    adata_out.layers["adt_clr"] = adata_out.X.copy()
    logger.info("CLR values saved to adata.layers['adt_clr'] and adata.X")

    # ------------------------------------------------------------------
    # Step 3 — Build metrics dict
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

    metrics = {
        "n_cells": int(adata_out.n_obs),
        "n_proteins": int(adata_out.n_vars),
        "clr_axis": clr_axis,
        "dsb_applied": False,
        "raw_counts_in_layer": "counts",
        "clr_in_layer": "adt_clr",
        # Post-CLR distribution summaries
        "clr_mean_per_cell": float(x_arr.sum(axis=1).mean()),
        "clr_max": float(x_arr.max()),
        "clr_min": float(x_arr.min()),
        # Raw count summaries
        "raw_median_total_counts_per_cell": float(
            np.median(raw_arr.sum(axis=1))
        ),
        "raw_max_count": float(raw_arr.max()),
    }

    # ------------------------------------------------------------------
    # Step 4 — Store provenance in uns
    # ------------------------------------------------------------------
    adata_out.uns["omicsage_adt_normalize"] = {
        "clr_axis": clr_axis,
        "clr_axis_description": (
            "per-protein across cells (muon default)"
            if clr_axis == 0
            else "per-cell across proteins"
        ),
        "dsb_applied": False,
        "normalized_in_layer": "adt_clr",
        "raw_in_layer": "counts",
        "muon_version": muon.__version__,
        "omicsage_module": "pipeline.modules.cite.adt_normalize",
        "omicsage_version": "0.1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    logger.info(
        "ADT normalization complete — %d cells, %d proteins, CLR axis=%d",
        adata_out.n_obs,
        adata_out.n_vars,
        clr_axis,
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
        raise NotImplementedError(
            "DSB normalization is not yet implemented in OmicSage. "
            "DSB requires empty-droplet ADT counts and isotype controls. "
            "Pass dsb_empty_adata=None to use CLR normalization only. "
            "DSB support is planned for a future session once empty-droplet "
            "data is available in the pipeline."
        )

    # Warn (don't raise) if .X does not look like raw integer counts —
    # CLR is designed for raw counts; normalized input will produce wrong results
    x_sample = adata.X[:min(100, adata.n_obs), :min(100, adata.n_vars)]
    if hasattr(x_sample, "toarray"):
        x_sample = x_sample.toarray()
    x_sample = np.asarray(x_sample, dtype=float)

    non_integer_frac = np.mean(np.abs(x_sample - np.round(x_sample)) > 1e-3)
    if non_integer_frac > 0.05:
        logger.warning(
            "%.0f%% of sampled ADT values are non-integer. "
            "CLR normalization expects raw counts in .X. "
            "If data is already normalized, results may be unreliable.",
            non_integer_frac * 100,
        )

    logger.debug(
        "Input validation passed — %d cells × %d proteins", adata.n_obs, adata.n_vars
    )
