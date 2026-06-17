"""
cite_deg.py — CITE-seq Differential Expression module for OmicSage
Phase 4, Step 7

Runs two complementary DEG analyses on CITE-seq data:

  A. DPE (Differential Protein Expression)
     Wilcoxon rank-sum test on the ADT layer, grouped by adt_celltype_manual
     (or adt_celltype_score as fallback). Identifies which proteins are
     significantly enriched per cell type — statistically confirming what
     the dotplot shows visually.

  B. Cross-modal RNA DEG
     Uses the SAME cell type grouping (ADT-defined) but runs DEG on the
     RNA layer (mdata["rna"].X / layers["logcounts"]). This finds the
     transcriptional programs underlying each immunophenotypically defined
     population. The gene lists produced here are fed directly into
     cite_gsea.py (Step 8).

Input
-----
Primary  : cite_06_integration.h5mu  (MuData — both RNA + ADT aligned)
Fallback : cite_05_annotated_adt.h5ad (AnnData — ADT only, DPE only)

When only the AnnData is available, cross-modal RNA DEG is skipped with a
UserWarning and gsea_cite (Step 8) will have nothing to process.

Outputs written to mdata (MuData path) or adata (AnnData fallback)
--------------------------------------------------------------------
uns["omicsage_cite_deg"]          Provenance block
uns["rank_genes_groups_dpe"]      Raw scanpy DPE results (ADT)
uns["rank_genes_groups_rna_cm"]   Raw scanpy cross-modal results (RNA)
                                  (only when MuData input is available)

Returns
-------
(MuData | AnnData, cite_deg_dict)

cite_deg_dict keys
------------------
  "dpe"          : {cell_type: DataFrame(protein, score, pval, logfc, pval_adj)}
  "rna_crossmodal: {cell_type: DataFrame(gene, score, pval, logfc, pval_adj)}
                   (empty dict when only AnnData fallback was used)
  "dpe_summary"  : wide DataFrame — top 5 proteins per cell type
  "rna_summary"  : wide DataFrame — top 5 genes per cell type
  "provenance"   : same metadata as uns["omicsage_cite_deg"]
  "input_type"   : "mudata" | "anndata"

Usage
-----
    from pipeline.modules.cite.cite_deg import cite_deg

    # MuData path (preferred — enables cross-modal DEG)
    mdata_deg, deg_dict = cite_deg(mdata, groupby="adt_celltype_manual")

    # AnnData fallback (DPE only)
    adata_deg, deg_dict = cite_deg(adata, groupby="adt_celltype_manual")
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

try:
    from mudata import MuData
    _MUDATA_AVAILABLE = True
except ImportError:
    _MUDATA_AVAILABLE = False
    MuData = None  # type: ignore


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cite_deg(
    data: Union["MuData", AnnData],
    groupby: str = "adt_celltype_manual",
    groupby_fallback: str = "adt_celltype_score",
    leiden_fallback: str = "leiden",
    method: str = "wilcoxon",
    min_logfc: float = 0.25,
    max_pval_adj: float = 0.05,
    n_genes: int = 200,
    exclude_protein_prefixes: Optional[list[str]] = None,
    exclude_gene_prefixes: Optional[list[str]] = None,
    use_raw_rna: bool = False,
    inplace: bool = False,
) -> tuple[Union["MuData", AnnData], dict]:
    """
    Run DPE (ADT layer) and cross-modal RNA DEG (RNA layer) on CITE-seq data.

    Parameters
    ----------
    data : MuData or AnnData
        MuData from cite_06 (preferred) or AnnData from cite_05 (fallback).
        MuData enables cross-modal RNA DEG; AnnData gives DPE only.
    groupby : str
        obs column on the ADT modality to group cells by.
        Default: "adt_celltype_manual" (from annotate_adt step).
    groupby_fallback : str
        Secondary fallback obs column if groupby is missing.
        Default: "adt_celltype_score".
    leiden_fallback : str
        Final fallback to raw Leiden cluster IDs if both above are missing.
        Default: "leiden".
    method : str
        Statistical test for rank_genes_groups. Default: "wilcoxon".
    min_logfc : float
        Minimum absolute log2 fold-change for a result to be retained.
    max_pval_adj : float
        Maximum BH-adjusted p-value for a result to be retained.
    n_genes : int
        Top N features to compute per group. Default: 200.
        For ADT with 134 proteins this effectively means all proteins.
    exclude_protein_prefixes : list of str, optional
        Protein name prefixes to exclude from DPE results after thresholding.
        Example: ["Mouse-IgG", "Rat-IgG"] to remove isotype controls.
        Default: None.
    exclude_gene_prefixes : list of str, optional
        Gene name prefixes to exclude from cross-modal RNA DEG results.
        Example: ["RPL", "RPS", "MT-"].
        Default: None.
    use_raw_rna : bool
        If False (default), uses mdata["rna"].layers["logcounts"] for RNA DEG.
        If True, uses mdata["rna"].X directly.
    inplace : bool
        If False (default), operates on a copy.

    Returns
    -------
    data : MuData or AnnData
        Input object with uns["omicsage_cite_deg"] populated.
    cite_deg_dict : dict
        See module docstring for key descriptions.

    Raises
    ------
    ValueError
        If no valid groupby column is found on the ADT modality.
    """
    is_mudata = _MUDATA_AVAILABLE and isinstance(data, MuData)

    if not inplace:
        data = data.copy()

    # Normalise prefix lists
    exclude_protein_prefixes = [
        p.upper() for p in (exclude_protein_prefixes or [])
    ]
    exclude_gene_prefixes = [
        p.upper() for p in (exclude_gene_prefixes or [])
    ]

    # ------------------------------------------------------------------
    # 1. Extract ADT AnnData and resolve groupby column
    # ------------------------------------------------------------------
    adt = data["adt"] if is_mudata else data

    group_col = _resolve_groupby(
        adt, groupby, groupby_fallback, leiden_fallback
    )

    # ------------------------------------------------------------------
    # 2. DPE — differential protein expression on ADT layer
    # ------------------------------------------------------------------
    dpe_results, dpe_summary = _run_dpe(
        adt=adt,
        group_col=group_col,
        method=method,
        min_logfc=min_logfc,
        max_pval_adj=max_pval_adj,
        n_genes=n_genes,
        exclude_prefixes=exclude_protein_prefixes,
    )

    # ------------------------------------------------------------------
    # 3. Cross-modal RNA DEG — RNA layer, same cell grouping
    # ------------------------------------------------------------------
    rna_results: dict = {}
    rna_summary: pd.DataFrame = pd.DataFrame()
    input_type = "anndata"

    if is_mudata:
        input_type = "mudata"
        rna_results, rna_summary = _run_crossmodal_rna_deg(
            mdata=data,
            group_col=group_col,
            method=method,
            min_logfc=min_logfc,
            max_pval_adj=max_pval_adj,
            n_genes=n_genes,
            exclude_prefixes=exclude_gene_prefixes,
            use_raw=use_raw_rna,
        )
    else:
        warnings.warn(
            "cite_deg received an AnnData (cite_05 fallback). "
            "Cross-modal RNA DEG requires MuData from cite_06. "
            "Only DPE will be computed. GSEA (cite_08) will be skipped.",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 4. Provenance
    # ------------------------------------------------------------------
    n_dpe_sig = sum(len(df) for df in dpe_results.values())
    n_rna_sig = sum(len(df) for df in rna_results.values())

    provenance = {
        "module": "cite_deg",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_type": input_type,
        "groupby": group_col,
        "method": method,
        "min_logfc": min_logfc,
        "max_pval_adj": max_pval_adj,
        "n_genes": n_genes,
        "exclude_protein_prefixes": exclude_protein_prefixes,
        "exclude_gene_prefixes": exclude_gene_prefixes,
        "n_cells": data.n_obs,
        "n_cell_types": len(dpe_results),
        "n_dpe_significant": n_dpe_sig,
        "n_rna_crossmodal_significant": n_rna_sig,
        "scanpy_version": sc.__version__,
        "outputs": {
            "dpe_key": "rank_genes_groups_dpe",
            "rna_cm_key": "rank_genes_groups_rna_cm" if is_mudata else None,
        },
    }

    # Write provenance to the top-level object
    if is_mudata:
        data.uns["omicsage_cite_deg"] = provenance
    else:
        data.uns["omicsage_cite_deg"] = provenance

    # ------------------------------------------------------------------
    # 5. Assemble return dict
    # ------------------------------------------------------------------
    cite_deg_dict = {
        "dpe": dpe_results,
        "rna_crossmodal": rna_results,
        "dpe_summary": dpe_summary,
        "rna_summary": rna_summary,
        "provenance": provenance,
        "input_type": input_type,
    }

    return data, cite_deg_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _resolve_groupby(
    adt: AnnData,
    groupby: str,
    groupby_fallback: str,
    leiden_fallback: str,
) -> str:
    """
    Return the first obs column found in priority order:
      groupby → groupby_fallback → leiden_fallback
    Raises ValueError if none are present.
    """
    for col in [groupby, groupby_fallback, leiden_fallback]:
        if col in adt.obs.columns:
            if col != groupby:
                warnings.warn(
                    f"Column '{groupby}' not found in adt.obs. "
                    f"Using '{col}' instead.",
                    UserWarning,
                    stacklevel=3,
                )
            return col

    raise ValueError(
        f"None of '{groupby}', '{groupby_fallback}', '{leiden_fallback}' "
        f"found in adt.obs. Available columns: {list(adt.obs.columns)}"
    )


def _run_dpe(
    adt: AnnData,
    group_col: str,
    method: str,
    min_logfc: float,
    max_pval_adj: float,
    n_genes: int,
    exclude_prefixes: list[str],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Run rank_genes_groups on the ADT layer (CLR-normalised).

    Uses adt.layers["adt_clr"] as the expression matrix.
    Results stored in adt.uns["rank_genes_groups_dpe"].

    Returns
    -------
    results : dict — {cell_type: DataFrame(protein, score, pval, logfc, pval_adj)}
    summary : DataFrame — top 5 proteins per cell type (long format)
    """
    # Warn on small groups
    group_counts = adt.obs[group_col].value_counts()
    small = group_counts[group_counts < 10].index.tolist()
    if small:
        warnings.warn(
            f"DPE: groups with < 10 cells (results may be unreliable): {small}",
            UserWarning,
            stacklevel=3,
        )

    # Set X to CLR layer for rank_genes_groups
    original_X = adt.X.copy()
    if "adt_clr" in adt.layers:
        adt.X = adt.layers["adt_clr"].copy()
    else:
        warnings.warn(
            "adt.layers['adt_clr'] not found — using adt.X as-is for DPE. "
            "Run adt_normalize first.",
            UserWarning,
            stacklevel=3,
        )

    try:
        sc.tl.rank_genes_groups(
            adt,
            groupby=group_col,
            method=method,
            n_genes=n_genes,
            rankby_abs=True,
            corr_method="benjamini-hochberg",
            use_raw=False,
            key_added="rank_genes_groups_dpe",
        )
    finally:
        # Always restore original X
        adt.X = original_X

    groups = list(adt.uns["rank_genes_groups_dpe"]["names"].dtype.names)
    results: dict[str, pd.DataFrame] = {}

    for group in groups:
        df = _extract_results(adt, group, n_genes, key="rank_genes_groups_dpe")
        df = _apply_thresholds(df, min_logfc, max_pval_adj)
        df = _apply_prefix_exclusion(df, exclude_prefixes, group,
                                     feature_col="gene")
        # Rename "gene" → "protein" for clarity in DPE output
        df = df.rename(columns={"gene": "protein"})
        results[group] = df

    summary = _build_summary(results, n_top=5, feature_col="protein")
    return results, summary


def _run_crossmodal_rna_deg(
    mdata: "MuData",
    group_col: str,
    method: str,
    min_logfc: float,
    max_pval_adj: float,
    n_genes: int,
    exclude_prefixes: list[str],
    use_raw: bool,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Run rank_genes_groups on the RNA layer using ADT-defined cell type labels.

    Transfers obs[group_col] from mdata["adt"] to mdata["rna"] using the
    shared barcode index, then runs DEG on the RNA layer.

    Results stored in mdata["rna"].uns["rank_genes_groups_rna_cm"].

    Returns
    -------
    results : dict — {cell_type: DataFrame(gene, score, pval, logfc, pval_adj)}
    summary : DataFrame — top 5 genes per cell type (long format)
    """
    rna = mdata["rna"]
    adt = mdata["adt"]

    # Transfer ADT cell type labels to RNA by shared barcode
    shared = rna.obs_names.intersection(adt.obs_names)
    if len(shared) == 0:
        warnings.warn(
            "Cross-modal DEG: RNA and ADT share no barcodes. "
            "Skipping cross-modal DEG.",
            UserWarning,
            stacklevel=3,
        )
        return {}, pd.DataFrame()

    rna_sub = rna[shared].copy()

    # Transfer groupby column from ADT obs to RNA obs
    rna_sub.obs["_cite_deg_group"] = (
        adt.obs.loc[shared, group_col].astype(str).values
    )

    # Prepare expression layer
    if not use_raw:
        if "logcounts" in rna_sub.layers:
            rna_sub.X = rna_sub.layers["logcounts"].copy()
        else:
            warnings.warn(
                "Cross-modal DEG: rna.layers['logcounts'] not found. "
                "Using rna.X as-is. Results may not be log-normalised.",
                UserWarning,
                stacklevel=3,
            )

    # Warn on small groups
    group_counts = rna_sub.obs["_cite_deg_group"].value_counts()
    small = group_counts[group_counts < 10].index.tolist()
    if small:
        warnings.warn(
            f"Cross-modal DEG: groups with < 10 cells: {small}",
            UserWarning,
            stacklevel=3,
        )

    sc.tl.rank_genes_groups(
        rna_sub,
        groupby="_cite_deg_group",
        method=method,
        n_genes=n_genes,
        rankby_abs=True,
        corr_method="benjamini-hochberg",
        use_raw=False,
        key_added="rank_genes_groups_rna_cm",
    )

    # Copy results key back to main RNA object
    mdata["rna"].uns["rank_genes_groups_rna_cm"] = (
        rna_sub.uns["rank_genes_groups_rna_cm"]
    )

    groups = list(rna_sub.uns["rank_genes_groups_rna_cm"]["names"].dtype.names)
    results: dict[str, pd.DataFrame] = {}

    for group in groups:
        df = _extract_results(
            rna_sub, group, n_genes, key="rank_genes_groups_rna_cm"
        )
        df = _apply_thresholds(df, min_logfc, max_pval_adj)
        df = _apply_prefix_exclusion(df, exclude_prefixes, group,
                                     feature_col="gene")
        results[group] = df

    summary = _build_summary(results, n_top=5, feature_col="gene")
    return results, summary


def _extract_results(
    adata: AnnData,
    group: str,
    n_genes: int,
    key: str = "rank_genes_groups",
) -> pd.DataFrame:
    """Extract rank_genes_groups results for one group into a tidy DataFrame."""
    rgg = adata.uns[key]

    def _safe(field: str) -> np.ndarray:
        if field not in rgg:
            return np.zeros(n_genes) if field != "names" else np.array([""] * n_genes)
        arr = rgg[field]
        if arr.dtype.names and group in arr.dtype.names:
            return arr[group]
        return np.zeros(n_genes) if field != "names" else np.array([""] * n_genes)

    df = pd.DataFrame({
        "gene":     _safe("names"),
        "score":    _safe("scores"),
        "pval":     _safe("pvals"),
        "logfc":    _safe("logfoldchanges"),
        "pval_adj": _safe("pvals_adj"),
    })
    df = df.dropna(subset=["gene"])
    df = df[df["gene"] != ""]
    return df.reset_index(drop=True)


def _apply_thresholds(
    df: pd.DataFrame,
    min_logfc: float,
    max_pval_adj: float,
) -> pd.DataFrame:
    mask = (df["pval_adj"] <= max_pval_adj) & (df["logfc"].abs() >= min_logfc)
    return df.loc[mask].sort_values("pval_adj").reset_index(drop=True)


def _apply_prefix_exclusion(
    df: pd.DataFrame,
    exclude_prefixes: list[str],
    group: str,
    feature_col: str = "gene",
) -> pd.DataFrame:
    if not exclude_prefixes or df.empty:
        return df
    mask_keep = ~df[feature_col].str.upper().apply(
        lambda g: any(g.startswith(p) for p in exclude_prefixes)
    )
    n_removed = (~mask_keep).sum()
    if n_removed > 0:
        removed = df.loc[~mask_keep, feature_col].tolist()
        warnings.warn(
            f"Group '{group}': excluded {n_removed} feature(s) matching "
            f"prefixes {exclude_prefixes}: {removed[:10]}"
            f"{'...' if len(removed) > 10 else ''}",
            UserWarning,
            stacklevel=4,
        )
    return df.loc[mask_keep].reset_index(drop=True)


def _build_summary(
    results: dict[str, pd.DataFrame],
    n_top: int = 5,
    feature_col: str = "gene",
) -> pd.DataFrame:
    """Build a long-format summary DataFrame with top N features per group."""
    rows = []
    for group, df in results.items():
        for rank, (_, row) in enumerate(df.head(n_top).iterrows(), start=1):
            rows.append({
                "group":    group,
                "rank":     rank,
                feature_col: row[feature_col],
                "logfc":    round(float(row["logfc"]), 4),
                "pval_adj": float(row["pval_adj"]),
            })
    return pd.DataFrame(rows)
