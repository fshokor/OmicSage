"""
adt_doublets.py — ADT doublet detection for OmicSage Phase 4

Detects heterotypic doublets using mutually exclusive cell-type surface
protein markers.  A real single cell cannot simultaneously express both a
T-cell marker (CD3) and a B-cell marker (CD19), so co-expression above a
CLR-value threshold is treated as evidence of a doublet.

Reference:
  https://www.sc-best-practices.org/surface_protein/doublet_detection.html

Input
-----
mdata["adt"].layers["adt_clr"]   — CLR-normalised ADT values (float)
                                   produced by adt_normalize.py

Outputs added to mdata["adt"].obs
---------------------------------
adt_doublet_score       float  fraction of marker pairs where both markers
                                exceed the threshold (0.0–1.0)
adt_predicted_doublet   bool   True when adt_doublet_score > 0

API
---
detect_adt_doublets(
    mdata,
    marker_pairs=None,
    threshold=2.5,
    filter_doublets=False,
    inplace=False,
)
→ (MuData, dict)

plot_doublet_scatter(mdata["adt"], pair=(marker_a, marker_b))
→ matplotlib Figure   (used by cite_doublets_report.py)
"""

from __future__ import annotations

import warnings
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Default marker pairs — mutually exclusive lineage markers
# (panel-agnostic: checked against available vars at runtime)
# ---------------------------------------------------------------------------
DEFAULT_MARKER_PAIRS: List[Tuple[str, str]] = [
    ("CD3", "CD19"),    # T cell vs B cell
    ("CD3", "CD14"),    # T cell vs Monocyte
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def detect_adt_doublets(
    mdata,
    marker_pairs: Optional[List[Tuple[str, str]]] = None,
    threshold: float = 2.5,
    filter_doublets: bool = False,
    inplace: bool = False,
):
    """Detect (and optionally remove) ADT heterotypic doublets.

    Parameters
    ----------
    mdata : MuData
        Must contain modality ``"adt"`` with ``layers["adt_clr"]``.
    marker_pairs : list of (str, str), optional
        Pairs of mutually exclusive markers.  Each string is matched
        against ``mdata["adt"].var_names`` by prefix (case-insensitive),
        so ``"CD3"`` matches ``"CD3-1"``, ``"CD3_TotalSeqB"``, etc.
        Defaults to ``[("CD3", "CD19"), ("CD3", "CD14")]``.
    threshold : float, default 2.5
        CLR expression level above which a marker is considered
        "expressed" for doublet calling.  2.5 is recommended by the
        sc-best-practices reference for CLR-normalised ADT data.
    filter_doublets : bool, default False
        If True, remove doublet-flagged cells from **both** ``mdata["adt"]``
        and ``mdata["rna"]`` (cross-modal sync).  When False, flags are
        written but no cells are removed.
    inplace : bool, default False
        If True, modify ``mdata`` in place and return ``(mdata, metrics)``.
        If False, work on a copy.

    Returns
    -------
    mdata : MuData
        Updated object (copy unless ``inplace=True``).
    metrics : dict
        Summary statistics — see Notes.

    Notes
    -----
    Metrics dict keys:
      n_cells_before        int   cells before any filtering
      n_doublets_detected   int   cells flagged as doublets
      pct_doublets          float percentage flagged
      pairs_evaluated       list  marker pairs actually evaluated
                                  (pairs with missing vars are skipped)
      pairs_skipped         list  pairs skipped due to missing vars
      n_cells_after         int   cells remaining (== n_cells_before when
                                  filter_doublets=False)
      threshold             float threshold used
      filter_doublets       bool  whether filtering was applied
      resolved_pairs        list  (idx_a, idx_b, name_a, name_b) for plotting
    """
    if marker_pairs is None:
        marker_pairs = DEFAULT_MARKER_PAIRS

    if not inplace:
        import copy
        new_mod = {"adt": mdata["adt"].copy()}
        if "rna" in mdata.mod:
            new_mod["rna"] = mdata["rna"].copy()
        mdata = copy.copy(mdata)
        mdata.mod = new_mod

    adata_adt = mdata["adt"]

    if "adt_clr" not in adata_adt.layers:
        raise KeyError(
            "mdata['adt'].layers['adt_clr'] not found. "
            "Run normalize_adt() before detect_adt_doublets()."
        )

    n_cells_before = adata_adt.n_obs
    clr = adata_adt.layers["adt_clr"]

    var_names = list(adata_adt.var_names)

    # ------------------------------------------------------------------
    # Resolve marker names → column indices
    # ------------------------------------------------------------------
    def _resolve_marker(name: str) -> Optional[int]:
        name_up = name.upper()
        for idx, vn in enumerate(var_names):
            if vn.upper().startswith(name_up):
                return idx
        return None

    pairs_evaluated: List[Tuple[str, str]] = []
    pairs_skipped: List[Tuple[str, str]] = []
    resolved: List[Tuple[int, int, str, str]] = []

    for (a, b) in marker_pairs:
        idx_a = _resolve_marker(a)
        idx_b = _resolve_marker(b)
        if idx_a is None or idx_b is None:
            pairs_skipped.append((a, b))
        else:
            resolved.append((idx_a, idx_b, a, b))
            pairs_evaluated.append((a, b))

    if not resolved:
        adata_adt.obs["adt_doublet_score"] = 0.0
        adata_adt.obs["adt_predicted_doublet"] = False
        metrics = _build_metrics(
            n_cells_before=n_cells_before,
            n_doublets=0,
            pairs_evaluated=[],
            pairs_skipped=pairs_skipped,
            n_cells_after=n_cells_before,
            threshold=threshold,
            filter_doublets=filter_doublets,
            resolved_pairs=[],
        )
        _write_provenance(adata_adt, metrics)
        return mdata, metrics

    # ------------------------------------------------------------------
    # Score each cell: fraction of pairs where both markers > threshold
    # ------------------------------------------------------------------
    n_pairs = len(resolved)
    pair_flags = np.zeros((n_cells_before, n_pairs), dtype=bool)

    for col_idx, (idx_a, idx_b, _a, _b) in enumerate(resolved):
        expr_a = _get_column(clr, idx_a)
        expr_b = _get_column(clr, idx_b)
        pair_flags[:, col_idx] = (expr_a > threshold) & (expr_b > threshold)

    doublet_score = pair_flags.sum(axis=1) / n_pairs
    predicted_doublet = doublet_score > 0

    adata_adt.obs["adt_doublet_score"] = doublet_score.astype(np.float64)
    adata_adt.obs["adt_predicted_doublet"] = predicted_doublet

    n_doublets = int(predicted_doublet.sum())

    # ------------------------------------------------------------------
    # Optional cross-modal filtering
    # ------------------------------------------------------------------
    n_cells_after = n_cells_before
    if filter_doublets and n_doublets > 0:
        keep_mask = ~predicted_doublet
        keep_barcodes = adata_adt.obs_names[keep_mask]

        mdata.mod["adt"] = adata_adt[keep_mask].copy()

        if "rna" in mdata.mod:
            rna = mdata["rna"]
            shared = rna.obs_names.isin(keep_barcodes)
            mdata.mod["rna"] = rna[shared].copy()

        n_cells_after = int(keep_mask.sum())

    metrics = _build_metrics(
        n_cells_before=n_cells_before,
        n_doublets=n_doublets,
        pairs_evaluated=pairs_evaluated,
        pairs_skipped=pairs_skipped,
        n_cells_after=n_cells_after,
        threshold=threshold,
        filter_doublets=filter_doublets,
        resolved_pairs=resolved,
    )

    _write_provenance(mdata["adt"], metrics)
    return mdata, metrics


# ---------------------------------------------------------------------------
# Tutorial scatter plot  (sc-best-practices ch. doublet_detection)
# ---------------------------------------------------------------------------

def plot_doublet_scatter(
    adt,
    pair: Optional[Tuple[str, str]] = None,
    threshold: float = 2.5,
    figsize: Tuple[float, float] = (6, 5),
):
    """
    Scatter plot of CLR expression for one mutually-exclusive marker pair,
    coloured by ``obs["adt_predicted_doublet"]``.

    Reproduces the sc-best-practices doublet-detection tutorial figure:
    each dot is a cell; axes are CLR expression of the two markers;
    the threshold lines mark the co-expression region; doublets are red.

    Parameters
    ----------
    adt : AnnData
        ``mdata["adt"]`` after ``detect_adt_doublets`` has been called.
        Must have ``layers["adt_clr"]`` and ``obs["adt_predicted_doublet"]``.
    pair : (str, str), optional
        Marker pair to plot (prefix match, case-insensitive).
        Defaults to the first evaluable pair from DEFAULT_MARKER_PAIRS.
    threshold : float
        CLR threshold line drawn on both axes.  Default: 2.5.
    figsize : (float, float)
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    if "adt_clr" not in adt.layers:
        raise KeyError("adt.layers['adt_clr'] not found.")

    if "adt_predicted_doublet" not in adt.obs.columns:
        raise KeyError(
            "adt.obs['adt_predicted_doublet'] not found. "
            "Run detect_adt_doublets() first."
        )

    # Resolve which pair to plot
    var_names = list(adt.var_names)

    def _resolve(name):
        name_up = name.upper()
        for idx, vn in enumerate(var_names):
            if vn.upper().startswith(name_up):
                return idx, vn
        return None, None

    if pair is None:
        pair = DEFAULT_MARKER_PAIRS[0]

    idx_a, vn_a = _resolve(pair[0])
    idx_b, vn_b = _resolve(pair[1])

    if idx_a is None or idx_b is None:
        # Fallback: try any pair from defaults
        for fallback in DEFAULT_MARKER_PAIRS[1:]:
            idx_a, vn_a = _resolve(fallback[0])
            idx_b, vn_b = _resolve(fallback[1])
            if idx_a is not None and idx_b is not None:
                break
        if idx_a is None:
            # Nothing resolvable — placeholder
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5,
                    f"Markers '{pair[0]}' / '{pair[1]}' not found in var_names.\n"
                    "Check marker names against the protein panel.",
                    ha="center", va="center", fontsize=10, color="#888",
                    transform=ax.transAxes, wrap=True)
            ax.axis("off")
            return fig

    clr = adt.layers["adt_clr"]
    expr_a = _get_column(clr, idx_a)
    expr_b = _get_column(clr, idx_b)
    doublet_flag = adt.obs["adt_predicted_doublet"].values.astype(bool)

    n_doublets  = int(doublet_flag.sum())
    n_singlets  = int((~doublet_flag).sum())
    n_total     = len(doublet_flag)
    pct         = 100.0 * n_doublets / max(n_total, 1)

    fig, ax = plt.subplots(figsize=figsize)

    # Singlets first (grey, behind)
    ax.scatter(
        expr_a[~doublet_flag], expr_b[~doublet_flag],
        s=3, alpha=0.4, color="#9ecae1", rasterized=True,
        label=f"Singlet (n={n_singlets:,})",
    )
    # Doublets on top (red)
    if n_doublets > 0:
        ax.scatter(
            expr_a[doublet_flag], expr_b[doublet_flag],
            s=5, alpha=0.75, color="#e63946", rasterized=True,
            label=f"Doublet (n={n_doublets:,}, {pct:.1f}%)",
        )

    # Threshold lines — define the co-expression quadrant
    ax.axvline(threshold, color="#e63946", linewidth=1.2, linestyle="--", alpha=0.8)
    ax.axhline(threshold, color="#e63946", linewidth=1.2, linestyle="--", alpha=0.8)

    # Shade the co-expression quadrant (both markers > threshold)
    xmax = max(expr_a.max(), threshold + 1)
    ymax = max(expr_b.max(), threshold + 1)
    ax.fill_between(
        [threshold, xmax], threshold, ymax,
        color="#e63946", alpha=0.06, label="Co-expression region",
    )

    ax.set_xlabel(f"CLR({vn_a})", fontsize=11)
    ax.set_ylabel(f"CLR({vn_b})", fontsize=11)
    ax.set_title(
        f"ADT Doublet Detection — {vn_a} vs {vn_b}\n"
        f"(threshold = {threshold}  |  {n_doublets:,} / {n_total:,} cells flagged)",
        fontsize=12, fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=9, markerscale=3)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_column(matrix, col_idx: int) -> np.ndarray:
    col = matrix[:, col_idx]
    if hasattr(col, "toarray"):
        return np.asarray(col.toarray()).ravel()
    return np.asarray(col).ravel()


def _build_metrics(
    *,
    n_cells_before: int,
    n_doublets: int,
    pairs_evaluated: list,
    pairs_skipped: list,
    n_cells_after: int,
    threshold: float,
    filter_doublets: bool,
    resolved_pairs: list,
) -> dict:
    pct = round(100.0 * n_doublets / n_cells_before, 3) if n_cells_before else 0.0
    return {
        "n_cells_before": n_cells_before,
        "n_doublets_detected": n_doublets,
        "pct_doublets": pct,
        "pairs_evaluated": [list(p) for p in pairs_evaluated],
        "pairs_skipped": [list(p) for p in pairs_skipped],
        "n_cells_after": n_cells_after,
        "threshold": threshold,
        "filter_doublets": filter_doublets,
        # Kept for report plotting — not serialised to h5ad provenance
        "_resolved_pairs": resolved_pairs,
    }


def _write_provenance(adata, metrics: dict) -> None:
    import datetime
    # Strip non-serialisable _resolved_pairs before writing to uns
    safe_metrics = {k: v for k, v in metrics.items() if not k.startswith("_")}
    adata.uns["omicsage_adt_doublets"] = {
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "omicsage_module": "adt_doublets",
        **safe_metrics,
    }
