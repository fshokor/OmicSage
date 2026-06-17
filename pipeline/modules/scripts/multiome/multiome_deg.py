"""
multiome_deg.py — RNA DEG and Differential Chromatin Accessibility for OmicSage.
Phase 5, Session 5

Runs two complementary differential analysis steps on multiome data:

  A. RNA DEG (Differential Gene Expression)
     Wilcoxon rank-sum test on the RNA layer, grouped by atac_celltype labels
     (transferred from RNA in atac_annotate.py).  Identifies transcriptional
     programs underlying each chromatin-accessibility-defined population.
     Gene lists produced here feed into downstream GSEA.

  B. DCA (Differential Chromatin Accessibility)
     Wilcoxon rank-sum test on the ATAC peak counts, same grouping.
     Identifies which genomic regions are differentially open per cell type.
     Results are in terms of peak IDs (chr:start-end).

Input
-----
Primary: MuData from multiome_integration.py (preferred — has joint embedding)
Fallback: MuData from atac_annotate.py (works fine — integration not required)

When only an AnnData is provided (ATAC only), RNA DEG is skipped with a
UserWarning and DCA runs alone.

Outputs written
---------------
mdata.uns["omicsage_multiome_deg"]              Provenance block
mdata["rna"].uns["rank_genes_groups_rna"]       Raw scanpy RNA DEG results
mdata["atac"].uns["rank_genes_groups_dca"]      Raw scanpy DCA results

Returns
-------
(MuData | AnnData, multiome_deg_dict)

multiome_deg_dict keys
----------------------
  "rna_deg"      : {cell_type: DataFrame(gene, score, pval, logfc, pval_adj)}
                   (empty dict when AnnData fallback or no RNA)
  "dca"          : {cell_type: DataFrame(peak, score, pval, logfc, pval_adj)}
  "rna_summary"  : DataFrame — top 5 genes per cell type (long format)
  "dca_summary"  : DataFrame — top 5 peaks per cell type (long format)
  "provenance"   : same metadata as uns["omicsage_multiome_deg"]
  "input_type"   : "mudata" | "anndata"

Usage
-----
    from pipeline.modules.multiome.multiome_deg import multiome_deg

    mdata_deg, deg_dict = multiome_deg(mdata, groupby="atac_celltype")

References
----------
sc-best-practices ATAC QC:
  https://www.sc-best-practices.org/chromatin_accessibility/quality_control.html
sc-best-practices paired integration:
  https://www.sc-best-practices.org/multimodal_integration/paired_integration.html
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

def multiome_deg(
    data: Union["MuData", AnnData],
    groupby: str = "atac_celltype",
    leiden_fallback: str = "atac_leiden",
    method: str = "wilcoxon",
    min_logfc: float = 0.25,
    max_pval_adj: float = 0.05,
    n_genes: int = 200,
    exclude_gene_prefixes: Optional[list[str]] = None,
    exclude_peak_prefixes: Optional[list[str]] = None,
    use_raw_rna: bool = False,
    inplace: bool = False,
) -> tuple[Union["MuData", AnnData], dict]:
    """
    Run RNA DEG and Differential Chromatin Accessibility on multiome data.

    Parameters
    ----------
    data : MuData or AnnData
        MuData with "rna" and "atac" modalities (preferred), or bare ATAC
        AnnData (DCA only, RNA DEG skipped).
    groupby : str
        obs column on the ATAC modality to group cells by.
        Default: "atac_celltype" (from atac_annotate step).
    leiden_fallback : str
        Fallback obs column if groupby is absent.
        Default: "atac_leiden".
    method : str
        Statistical test for rank_genes_groups. Default: "wilcoxon".
    min_logfc : float
        Minimum absolute log2 fold-change threshold. Default: 0.25.
    max_pval_adj : float
        Maximum BH-adjusted p-value threshold. Default: 0.05.
    n_genes : int
        Top N features per group. Default: 200.
    exclude_gene_prefixes : list of str, optional
        Gene name prefixes to exclude from RNA DEG results.
        Example: ["RPL", "RPS", "MT-"].
    exclude_peak_prefixes : list of str, optional
        Peak name prefixes to exclude from DCA results.
        Example: ["chrM"] to drop mitochondrial peaks.
    use_raw_rna : bool
        If False (default), uses mdata["rna"].layers["logcounts"] for RNA DEG.
        If True, uses mdata["rna"].X directly.
    inplace : bool
        If False (default), operates on a copy.

    Returns
    -------
    data : MuData or AnnData
        Input object with provenance written to uns.
    multiome_deg_dict : dict
        See module docstring for key descriptions.

    Raises
    ------
    ValueError
        If no valid groupby column is found on the ATAC modality.
    """
    is_mudata = _MUDATA_AVAILABLE and isinstance(data, MuData)

    if not inplace:
        data = data.copy()

    exclude_gene_prefixes  = [p.upper() for p in (exclude_gene_prefixes or [])]
    exclude_peak_prefixes  = [p.upper() for p in (exclude_peak_prefixes or [])]

    # ------------------------------------------------------------------
    # 1. Extract ATAC AnnData and resolve groupby column
    # ------------------------------------------------------------------
    atac = data["atac"] if is_mudata else data

    group_col = _resolve_groupby(atac, groupby, leiden_fallback)

    # ------------------------------------------------------------------
    # 2. RNA DEG — RNA layer, ATAC-defined cell type grouping
    # ------------------------------------------------------------------
    rna_results: dict = {}
    rna_summary: pd.DataFrame = pd.DataFrame()
    input_type = "anndata"

    if is_mudata:
        input_type = "mudata"
        rna_results, rna_summary = _run_rna_deg(
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
            "multiome_deg received an AnnData (ATAC only). "
            "RNA DEG requires MuData with both 'rna' and 'atac' modalities. "
            "Only DCA will be computed.",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 3. DCA — differential chromatin accessibility on ATAC peaks
    # ------------------------------------------------------------------
    dca_results, dca_summary = _run_dca(
        atac=atac,
        group_col=group_col,
        method=method,
        min_logfc=min_logfc,
        max_pval_adj=max_pval_adj,
        n_genes=n_genes,
        exclude_prefixes=exclude_peak_prefixes,
    )

    # ------------------------------------------------------------------
    # 4. Provenance
    # ------------------------------------------------------------------
    n_rna_sig = sum(len(df) for df in rna_results.values())
    n_dca_sig = sum(len(df) for df in dca_results.values())

    provenance = {
        "module":    "multiome_deg",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_type": input_type,
        "groupby":   group_col,
        "method":    method,
        "min_logfc": min_logfc,
        "max_pval_adj": max_pval_adj,
        "n_genes":   n_genes,
        "exclude_gene_prefixes":  exclude_gene_prefixes,
        "exclude_peak_prefixes":  exclude_peak_prefixes,
        "n_cells":       int(data.n_obs),
        "n_cell_types":  int(len(dca_results)),
        "n_rna_significant": int(n_rna_sig),
        "n_dca_significant": int(n_dca_sig),
        "scanpy_version": sc.__version__,
        "outputs": {
            "rna_deg_key": "rank_genes_groups_rna" if is_mudata else None,
            "dca_key":     "rank_genes_groups_dca",
        },
    }

    data.uns["omicsage_multiome_deg"] = provenance

    # ------------------------------------------------------------------
    # 5. Assemble return dict
    # ------------------------------------------------------------------
    multiome_deg_dict = {
        "rna_deg":     rna_results,
        "dca":         dca_results,
        "rna_summary": rna_summary,
        "dca_summary": dca_summary,
        "provenance":  provenance,
        "input_type":  input_type,
    }

    return data, multiome_deg_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _resolve_groupby(
    atac: AnnData,
    groupby: str,
    leiden_fallback: str,
) -> str:
    """
    Return the first obs column found in priority order:
      groupby → leiden_fallback
    Raises ValueError if neither is present.
    """
    for col in [groupby, leiden_fallback]:
        if col in atac.obs.columns:
            if col != groupby:
                warnings.warn(
                    f"Column '{groupby}' not found in atac.obs. "
                    f"Using '{col}' instead.",
                    UserWarning,
                    stacklevel=3,
                )
            return col

    raise ValueError(
        f"None of '{groupby}', '{leiden_fallback}' "
        f"found in atac.obs. Available columns: {list(atac.obs.columns)}"
    )


def _run_rna_deg(
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
    Run rank_genes_groups on the RNA layer using ATAC-defined cell type labels.

    Transfers obs[group_col] from mdata["atac"] to mdata["rna"] via shared
    barcodes, then runs DEG on the RNA layer.

    Results stored in mdata["rna"].uns["rank_genes_groups_rna"].

    Returns
    -------
    results : dict — {cell_type: DataFrame(gene, score, pval, logfc, pval_adj)}
    summary : DataFrame — top 5 genes per cell type (long format)
    """
    rna  = mdata["rna"]
    atac = mdata["atac"]

    shared = rna.obs_names.intersection(atac.obs_names)
    if len(shared) == 0:
        warnings.warn(
            "RNA DEG: RNA and ATAC share no barcodes. Skipping RNA DEG.",
            UserWarning,
            stacklevel=3,
        )
        return {}, pd.DataFrame()

    rna_sub = rna[shared].copy()

    # Transfer groupby column from ATAC obs to RNA obs
    rna_sub.obs["_multiome_deg_group"] = (
        atac.obs.loc[shared, group_col].astype(str).values
    )

    # Prepare expression layer
    if not use_raw:
        if "logcounts" in rna_sub.layers:
            rna_sub.X = rna_sub.layers["logcounts"].copy()
        else:
            warnings.warn(
                "RNA DEG: rna.layers['logcounts'] not found. "
                "Using rna.X as-is. Results may not be log-normalised.",
                UserWarning,
                stacklevel=3,
            )

    group_counts = rna_sub.obs["_multiome_deg_group"].value_counts()
    small = group_counts[group_counts < 10].index.tolist()
    if small:
        warnings.warn(
            f"RNA DEG: groups with < 10 cells (results may be unreliable): {small}",
            UserWarning,
            stacklevel=3,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            rna_sub,
            groupby="_multiome_deg_group",
            method=method,
            n_genes=n_genes,
            rankby_abs=True,
            corr_method="benjamini-hochberg",
            use_raw=False,
            key_added="rank_genes_groups_rna",
        )

    # Copy results key back to the main RNA object
    mdata["rna"].uns["rank_genes_groups_rna"] = (
        rna_sub.uns["rank_genes_groups_rna"]
    )

    groups = list(rna_sub.uns["rank_genes_groups_rna"]["names"].dtype.names)
    results: dict[str, pd.DataFrame] = {}

    for group in groups:
        df = _extract_results(rna_sub, group, n_genes, key="rank_genes_groups_rna")
        df = _apply_thresholds(df, min_logfc, max_pval_adj)
        df = _apply_prefix_exclusion(df, exclude_prefixes, group, feature_col="gene")
        results[group] = df

    summary = _build_summary(results, n_top=5, feature_col="gene")
    return results, summary


def _run_dca(
    atac: AnnData,
    group_col: str,
    method: str,
    min_logfc: float,
    max_pval_adj: float,
    n_genes: int,
    exclude_prefixes: list[str],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Run rank_genes_groups on ATAC peak counts for differential accessibility.

    Uses raw peak counts from atac.layers["counts"] when available, otherwise
    atac.X (TF-IDF). TF-IDF is acceptable for rank-based tests but raw counts
    are preferred for fold-change interpretation.

    Results stored in atac.uns["rank_genes_groups_dca"].

    Returns
    -------
    results : dict — {cell_type: DataFrame(peak, score, pval, logfc, pval_adj)}
    summary : DataFrame — top 5 peaks per cell type (long format)
    """
    group_counts = atac.obs[group_col].value_counts()
    small = group_counts[group_counts < 10].index.tolist()
    if small:
        warnings.warn(
            f"DCA: groups with < 10 cells (results may be unreliable): {small}",
            UserWarning,
            stacklevel=3,
        )

    # Swap in raw counts for DCA if available; preserve original X
    original_X = atac.X.copy() if hasattr(atac.X, "copy") else atac.X
    use_counts = "counts" in atac.layers

    if use_counts:
        atac.X = atac.layers["counts"].copy()
    else:
        warnings.warn(
            "atac.layers['counts'] not found — using atac.X (TF-IDF) for DCA. "
            "Fold-change values will be on the TF-IDF scale, not raw counts. "
            "For best results, ensure atac_qc.py preserved raw counts in layers['counts'].",
            UserWarning,
            stacklevel=3,
        )

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.rank_genes_groups(
                atac,
                groupby=group_col,
                method=method,
                n_genes=n_genes,
                rankby_abs=True,
                corr_method="benjamini-hochberg",
                use_raw=False,
                key_added="rank_genes_groups_dca",
            )
    finally:
        atac.X = original_X

    groups = list(atac.uns["rank_genes_groups_dca"]["names"].dtype.names)
    results: dict[str, pd.DataFrame] = {}

    for group in groups:
        df = _extract_results(atac, group, n_genes, key="rank_genes_groups_dca")
        df = _apply_thresholds(df, min_logfc, max_pval_adj)
        df = _apply_prefix_exclusion(df, exclude_prefixes, group, feature_col="gene")
        # Rename "gene" → "peak" for clarity in DCA output
        df = df.rename(columns={"gene": "peak"})
        results[group] = df

    summary = _build_summary(results, n_top=5, feature_col="peak")
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
            return (
                np.zeros(n_genes) if field != "names"
                else np.array([""] * n_genes)
            )
        arr = rgg[field]
        if arr.dtype.names and group in arr.dtype.names:
            return arr[group]
        return (
            np.zeros(n_genes) if field != "names"
            else np.array([""] * n_genes)
        )

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
    n_removed = int((~mask_keep).sum())
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
                "group":      group,
                "rank":       rank,
                feature_col:  row[feature_col],
                "logfc":      round(float(row["logfc"]), 4),
                "pval_adj":   float(row["pval_adj"]),
            })
    return pd.DataFrame(rows)
