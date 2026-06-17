"""
cite_epitope.py — Epitope Characterisation module for OmicSage
Phase 4, Step 10

Four sub-analyses that together replicate and extend the Seurat multimodal
vignette (https://satijalab.org/seurat/articles/multimodal_vignette):

  A. Epitope abundance scoring
     Per-cell lineage panel scores computed from CLR-normalised ADT values.
     Writes new obs columns: obs["epitope_score_<panel_name>"].
     Two panel sources (custom > preset, same pattern as adt_annotate.py):
       - preset="bmmc"         : built-in BioLegend TotalSeq-A panel
       - epitope_panels={...}  : user-supplied {panel_name: [protein, ...]}

  B. Ridge plots — ADT markers per cell type
     CLR distribution of the top marker proteins per ADT-defined cell type,
     as overlapping ridgeline plots coloured by cell type. Shows how cleanly
     each protein separates its target population.

  C. RNA-protein parallel plots
     For matched protein-gene pairs (from cite_09 correlation results),
     side-by-side violin plots: CLR expression on the left, log-normalised
     RNA on the right, same cell type grouping. Makes RNA-protein divergence
     visible per cell type — the core Seurat multimodal vignette figure.
     Requires corr_results DataFrame from cite_09. Skipped gracefully if absent.

  D. Cross-modal marker table
     Per cell type: top ADT markers (from DPE in cite_07) annotated with
     their RNA counterpart gene and Spearman r from cite_09.
     Synthesis output that ties cite_05 + cite_07 + cite_09 together.
     Requires dpe_results from cite_07. Skipped gracefully if absent.

Input
-----
  mdata          : MuData from cite_06 — required
  corr_results   : DataFrame from cite_09 corr_dict["results"] — optional
  dpe_results    : dict from cite_07 deg_dict["dpe"] — optional
  preset         : str — "bmmc" or any key in _PRESET_EPITOPE_PANELS
  epitope_panels : dict[str, list[str]] — custom panels, takes precedence

Outputs written to mdata
------------------------
  obs["epitope_score_<panel>"]   float  Per-cell score per panel (always)
  uns["omicsage_cite_epitope"]          Provenance block

Returns
-------
(MuData, epitope_dict)

epitope_dict keys
-----------------
  "scores_df"    : DataFrame(cell, panel_name, score) — long format
  "panel_used"   : dict[str, list[str]] — resolved panels
  "n_panels"     : int
  "provenance"   : same as uns block
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

try:
    from mudata import MuData
    _MUDATA_AVAILABLE = True
except ImportError:
    _MUDATA_AVAILABLE = False
    MuData = None  # type: ignore


# ---------------------------------------------------------------------------
# Built-in epitope panels
# Same panel names as adt_annotate._PRESET_PANELS["bmmc"] so the two modules
# are consistent — custom > preset, same resolution pattern.
# ---------------------------------------------------------------------------

_PRESET_EPITOPE_PANELS: dict[str, dict[str, list[str]]] = {
    "bmmc": {
        # Groups of proteins that define each lineage.
        # Protein names are prefix-matched against adt.var_names
        # (case-insensitive) so "CD3" matches "CD3E", "CD3D", etc.
        "T_cell":    ["CD3", "CD5", "CD2"],
        "CD4_T":     ["CD4", "CD3", "CD5"],
        "CD8_T":     ["CD8a", "CD3", "CD5"],
        "NK":        ["CD56", "CD16", "NKG7"],
        "B_cell":    ["CD19", "CD20", "CD22"],
        "Myeloid":   ["CD14", "CD11c", "CD172a"],
        "Erythroid": ["CD71", "CD36", "CD235a"],
        "HSC":       ["CD34", "CD38", "CD90"],
        "Plasma":    ["CD38", "CD138"],
    },
    "pbmc": {
        "T_cell":  ["CD3", "CD5"],
        "CD4_T":   ["CD4", "CD3"],
        "CD8_T":   ["CD8a", "CD3"],
        "NK":      ["CD56", "CD16"],
        "B_cell":  ["CD19", "CD20"],
        "Mono":    ["CD14", "CD16"],
        "DC":      ["CD123", "CD11c"],
    },
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cite_epitope(
    mdata: "MuData",
    corr_results: Optional[pd.DataFrame] = None,
    dpe_results: Optional[dict[str, pd.DataFrame]] = None,
    preset: Optional[str] = None,
    epitope_panels: Optional[dict[str, list[str]]] = None,
    groupby: str = "adt_celltype_manual",
    groupby_fallback: str = "adt_celltype_score",
    leiden_fallback: str = "leiden",
    n_top_markers: int = 5,
    inplace: bool = False,
) -> tuple["MuData", dict]:
    """
    Run epitope characterisation on CITE-seq data.

    Parameters
    ----------
    mdata : MuData
        MuData from cite_06. Requires mdata["adt"] and mdata["rna"].
    corr_results : DataFrame or None
        DataFrame from cite_09 corr_dict["results"].
        Columns: protein, gene, r, pval, pval_adj, n_cells, matched_by.
        When None, sub-analysis C (RNA-protein parallel) and the r column
        in sub-analysis D are skipped.
    dpe_results : dict or None
        Dict from cite_07 deg_dict["dpe"]:
        {cell_type: DataFrame(protein, score, pval, logfc, pval_adj)}.
        When None, sub-analysis D (cross-modal marker table) is skipped.
    preset : str or None
        Built-in panel name. Currently: "bmmc", "pbmc".
        Used when epitope_panels is None.
    epitope_panels : dict[str, list[str]] or None
        Custom panels: {panel_name: [protein_prefix, ...]}.
        Takes precedence over preset.
        Protein names are prefix-matched against adt.var_names
        (case-insensitive).
    groupby : str
        obs column on mdata["adt"] for cell type grouping.
        Default: "adt_celltype_manual".
    groupby_fallback : str
        Secondary fallback. Default: "adt_celltype_score".
    leiden_fallback : str
        Final fallback to raw Leiden cluster IDs.
    n_top_markers : int
        Number of top ADT markers per cell type used in the cross-modal
        marker table (sub-analysis D). Default: 5.
    inplace : bool
        If False (default), operates on a copy.

    Returns
    -------
    mdata : MuData
        Input MuData with obs["epitope_score_*"] columns and
        uns["omicsage_cite_epitope"] populated.
    epitope_dict : dict
        Keys: scores_df, panel_used, n_panels, marker_table, provenance.
    """
    if not (_MUDATA_AVAILABLE and isinstance(mdata, MuData)):
        raise ValueError(
            "cite_epitope requires a MuData object with 'rna' and 'adt' "
            "modalities (cite_06 checkpoint)."
        )

    if not inplace:
        mdata = mdata.copy()

    adt = mdata["adt"]

    # ------------------------------------------------------------------
    # 1. Resolve groupby column
    # ------------------------------------------------------------------
    group_col = _resolve_groupby(adt, groupby, groupby_fallback, leiden_fallback)

    # ------------------------------------------------------------------
    # 2. Resolve epitope panels  (custom > preset > warn + empty)
    # ------------------------------------------------------------------
    active_panels = _resolve_panels(epitope_panels, preset)

    # ------------------------------------------------------------------
    # 3. Sub-analysis A — Epitope abundance scoring
    # ------------------------------------------------------------------
    scores_df, score_keys = _compute_epitope_scores(adt, active_panels)

    # Write per-cell scores as obs columns
    for panel_name, col_key in score_keys.items():
        adt.obs[col_key] = scores_df.set_index("cell")[panel_name].values

    # ------------------------------------------------------------------
    # 4. Sub-analysis D — Cross-modal marker table
    #    (built from dpe_results + corr_results, both optional)
    # ------------------------------------------------------------------
    marker_table = _build_marker_table(
        dpe_results=dpe_results,
        corr_results=corr_results,
        n_top=n_top_markers,
    )

    # ------------------------------------------------------------------
    # 5. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "module":           "cite_epitope",
        "timestamp":        datetime.now(timezone.utc).isoformat(),
        "groupby":          group_col,
        "preset":           preset,
        "n_panels":         len(active_panels),
        "panel_names":      list(active_panels.keys()),
        "score_obs_keys":   list(score_keys.values()),
        "has_corr_results": corr_results is not None and not corr_results.empty,
        "has_dpe_results":  dpe_results is not None and len(dpe_results) > 0,
        "n_top_markers":    n_top_markers,
    }
    mdata.uns["omicsage_cite_epitope"] = provenance

    epitope_dict = {
        "scores_df":   scores_df,
        "panel_used":  active_panels,
        "n_panels":    len(active_panels),
        "marker_table": marker_table,
        "provenance":  provenance,
    }

    return mdata, epitope_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _resolve_groupby(
    adt: AnnData,
    groupby: str,
    groupby_fallback: str,
    leiden_fallback: str,
) -> str:
    for col in [groupby, groupby_fallback, leiden_fallback]:
        if col in adt.obs.columns:
            if col != groupby:
                warnings.warn(
                    f"Column '{groupby}' not found in adt.obs. "
                    f"Using '{col}' instead.",
                    UserWarning, stacklevel=3,
                )
            return col
    raise ValueError(
        f"None of '{groupby}', '{groupby_fallback}', '{leiden_fallback}' "
        f"found in adt.obs. Available: {list(adt.obs.columns)}"
    )


def _resolve_panels(
    epitope_panels: Optional[dict],
    preset: Optional[str],
) -> dict[str, list[str]]:
    """
    Return active panels. Resolution order: custom > preset > empty + warn.
    """
    if epitope_panels is not None:
        return {k: list(v) for k, v in epitope_panels.items()}

    if preset is not None:
        key = preset.lower().strip()
        if key in _PRESET_EPITOPE_PANELS:
            return dict(_PRESET_EPITOPE_PANELS[key])
        warnings.warn(
            f"Unknown preset '{preset}'. "
            f"Available: {list(_PRESET_EPITOPE_PANELS.keys())}. "
            "Epitope scoring will be skipped.",
            UserWarning, stacklevel=3,
        )

    warnings.warn(
        "No epitope_panels or preset provided. "
        "Epitope scoring will be skipped. "
        "Set preset='bmmc' or pass epitope_panels={{...}} in the config.",
        UserWarning, stacklevel=3,
    )
    return {}


def _resolve_proteins(
    panel_proteins: list[str],
    var_names: list[str],
) -> list[str]:
    """
    Resolve panel protein prefixes to actual var_names.
    Case-insensitive prefix match: "CD3" matches "CD3E", "CD3D", "CD3G".
    Returns list of matched var_names (may be empty if panel uses names
    not present in this dataset's antibody cocktail).
    """
    var_lower = [(v.lower(), v) for v in var_names]
    resolved: list[str] = []
    for prefix in panel_proteins:
        p_lower = prefix.lower()
        for vl, v_orig in var_lower:
            if vl.startswith(p_lower) and v_orig not in resolved:
                resolved.append(v_orig)
    return resolved


def _compute_epitope_scores(
    adt: AnnData,
    panels: dict[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Compute per-cell mean CLR score for each panel.

    Score = mean CLR expression of resolved panel proteins per cell.
    Cells where no panel proteins are detected get score = 0.

    Returns
    -------
    scores_df : DataFrame(cell, panel_name_1, panel_name_2, ...)
    score_keys : dict[panel_name, obs_col_name]
    """
    if not panels:
        empty_df = pd.DataFrame({"cell": adt.obs_names.tolist()})
        return empty_df, {}

    # Get CLR matrix
    if "adt_clr" in adt.layers:
        mat = adt.layers["adt_clr"]
    else:
        warnings.warn(
            "adt.layers['adt_clr'] not found — using adt.X for epitope scoring.",
            UserWarning, stacklevel=3,
        )
        mat = adt.X

    mat = mat.toarray() if sp.issparse(mat) else np.asarray(mat, dtype=float)
    var_names = list(adt.var_names)

    scores: dict[str, np.ndarray] = {}
    score_keys: dict[str, str] = {}

    for panel_name, protein_prefixes in panels.items():
        resolved = _resolve_proteins(protein_prefixes, var_names)

        if not resolved:
            warnings.warn(
                f"Panel '{panel_name}': no proteins matched var_names. "
                f"Requested prefixes: {protein_prefixes}. "
                "Score will be all zeros.",
                UserWarning, stacklevel=3,
            )
            scores[panel_name] = np.zeros(adt.n_obs)
        else:
            idxs = [var_names.index(p) for p in resolved]
            scores[panel_name] = mat[:, idxs].mean(axis=1)

        # obs column key: "epitope_score_T_cell" etc.
        safe_name = panel_name.replace(" ", "_").replace("-", "_")
        score_keys[panel_name] = f"epitope_score_{safe_name}"

    scores_df = pd.DataFrame({"cell": adt.obs_names.tolist(), **scores})
    return scores_df, score_keys


def _build_marker_table(
    dpe_results: Optional[dict[str, pd.DataFrame]],
    corr_results: Optional[pd.DataFrame],
    n_top: int,
) -> pd.DataFrame:
    """
    Build cross-modal marker table: top ADT markers per cell type,
    annotated with RNA counterpart and Spearman r.

    Returns empty DataFrame if dpe_results is None.
    """
    if not dpe_results:
        return pd.DataFrame(
            columns=["cell_type", "protein", "logfc", "pval_adj",
                     "rna_gene", "r", "r_pval_adj"]
        )

    # Build protein → r lookup from corr_results
    corr_lookup: dict[str, tuple[str, float, float]] = {}
    if corr_results is not None and not corr_results.empty:
        for _, row in corr_results.iterrows():
            corr_lookup[str(row["protein"])] = (
                str(row["gene"]),
                float(row["r"]),
                float(row.get("pval_adj", float("nan"))),
            )

    rows: list[dict] = []
    for cell_type, df in dpe_results.items():
        if df.empty:
            continue
        top = df.head(n_top)
        for _, row in top.iterrows():
            protein = str(row["protein"])
            gene, r, r_padj = corr_lookup.get(protein, ("—", float("nan"), float("nan")))
            rows.append({
                "cell_type": cell_type,
                "protein":   protein,
                "logfc":     float(row["logfc"]),
                "pval_adj":  float(row["pval_adj"]),
                "rna_gene":  gene,
                "r":         r,
                "r_pval_adj": r_padj,
            })

    return pd.DataFrame(rows)
