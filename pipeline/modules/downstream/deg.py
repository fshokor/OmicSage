"""
deg.py — Differential Expression Gene module for OmicSage
Phase 1, Step 6

Performs Wilcoxon rank-sum DEG analysis via scanpy.
Supports one-vs-rest (default) and pairwise modes.
Results written to adata.uns['omicsage_deg'] and returned as deg_dict.

Usage:
    from pipeline.modules.qc.deg import deg

    adata_deg, deg_dict = deg(
        adata_annotated,
        groupby="cell_type_vote",
        leiden_col="leiden",
        method="wilcoxon",
        min_logfc=0.25,
        max_pval_adj=0.05,
        n_genes=200,
        use_raw=False,
        inplace=False,
    )
"""

from __future__ import annotations

import warnings
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def deg(
    adata: AnnData,
    groupby: str = "cell_type_vote",
    leiden_col: str = "leiden",
    method: str = "wilcoxon",
    min_logfc: float = 0.25,
    max_pval_adj: float = 0.05,
    n_genes: int = 200,
    use_raw: bool = False,
    pairwise_groups: Optional[list[tuple[str, str]]] = None,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Run differential expression analysis on an annotated AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated AnnData. Should contain obs[groupby] or obs[leiden_col].
        Expression values must be in adata.layers['logcounts'] (use_raw=False)
        or adata.X / adata.raw (use_raw=True).
    groupby : str
        obs column to group cells by (default: 'cell_type_vote').
        Falls back to leiden_col if not present.
    leiden_col : str
        Fallback obs column if groupby is missing (default: 'leiden').
    method : str
        DEG method passed to sc.tl.rank_genes_groups.
        One of: 'wilcoxon', 't-test', 'logreg'.
    min_logfc : float
        Minimum absolute log2 fold-change threshold for filtering results.
    max_pval_adj : float
        Maximum adjusted p-value (BH FDR) threshold for filtering results.
    n_genes : int
        Number of top genes to compute per group.
    use_raw : bool
        If False (default), uses adata.layers['logcounts'].
        If True, uses adata.raw or adata.X.
    pairwise_groups : list of (str, str) tuples, optional
        If provided, also runs pairwise DEG for each specified pair.
        Each tuple is (group_a, group_b).
    inplace : bool
        If False (default), operates on a copy of adata.
        If True, modifies adata in place.

    Returns
    -------
    adata : AnnData
        AnnData with adata.uns['omicsage_deg'] populated.
    deg_dict : dict
        Dictionary with keys:
            'results'    — {group: DataFrame(gene, score, pval, logfc, pval_adj)}
            'summary_df' — wide DataFrame: top 5 genes per group
            'provenance' — same metadata as uns
            'pairwise'   — {(a, b): DataFrame} if pairwise_groups provided

    Raises
    ------
    ValueError
        If neither groupby nor leiden_col is found in adata.obs.
    """
    if not inplace:
        adata = adata.copy()

    # ------------------------------------------------------------------
    # 1. Resolve groupby column
    # ------------------------------------------------------------------
    if groupby in adata.obs.columns:
        group_col = groupby
    elif leiden_col in adata.obs.columns:
        warnings.warn(
            f"Column '{groupby}' not found in adata.obs. "
            f"Falling back to '{leiden_col}'.",
            UserWarning,
            stacklevel=2,
        )
        group_col = leiden_col
    else:
        raise ValueError(
            f"Neither '{groupby}' nor '{leiden_col}' found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # ------------------------------------------------------------------
    # 2. Warn about small groups (< 10 cells)
    # ------------------------------------------------------------------
    group_counts = adata.obs[group_col].value_counts()
    small_groups = group_counts[group_counts < 10].index.tolist()
    if small_groups:
        warnings.warn(
            f"The following groups have fewer than 10 cells and may produce "
            f"unreliable DEG results: {small_groups}",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 3. Prepare expression layer
    # ------------------------------------------------------------------
    if not use_raw:
        if "logcounts" not in adata.layers:
            raise ValueError(
                "use_raw=False requires adata.layers['logcounts']. "
                "Run normalize() before deg()."
            )
        # Temporarily set X to logcounts for rank_genes_groups
        adata.X = adata.layers["logcounts"]

    # ------------------------------------------------------------------
    # 4. Run one-vs-rest DEG
    # ------------------------------------------------------------------
    sc.tl.rank_genes_groups(
        adata,
        groupby=group_col,
        method=method,
        n_genes=n_genes,
        corr_method="benjamini-hochberg",
        use_raw=use_raw,
        pts=True,           # compute fraction of cells expressing gene
    )

    # ------------------------------------------------------------------
    # 5. Extract results into DataFrames
    # ------------------------------------------------------------------
    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    results = {}

    for group in groups:
        df = _extract_group_results(adata, group, n_genes)
        df = _apply_thresholds(df, min_logfc, max_pval_adj)
        results[group] = df

    # ------------------------------------------------------------------
    # 6. Optional pairwise DEG
    # ------------------------------------------------------------------
    pairwise_results = {}
    if pairwise_groups:
        for group_a, group_b in pairwise_groups:
            pair_key = (group_a, group_b)
            try:
                sc.tl.rank_genes_groups(
                    adata,
                    groupby=group_col,
                    groups=[group_a],
                    reference=group_b,
                    method=method,
                    n_genes=n_genes,
                    corr_method="benjamini-hochberg",
                    use_raw=use_raw,
                    key_added=f"rank_genes_pairwise_{group_a}_vs_{group_b}",
                )
                df = _extract_group_results(
                    adata,
                    group_a,
                    n_genes,
                    key=f"rank_genes_pairwise_{group_a}_vs_{group_b}",
                )
                df = _apply_thresholds(df, min_logfc, max_pval_adj)
                pairwise_results[pair_key] = df
            except Exception as e:
                warnings.warn(
                    f"Pairwise DEG failed for {group_a} vs {group_b}: {e}",
                    UserWarning,
                    stacklevel=2,
                )

    # ------------------------------------------------------------------
    # 7. Build summary DataFrame (top 5 genes per group)
    # ------------------------------------------------------------------
    summary_df = _build_summary_df(results, n_top=5)

    # ------------------------------------------------------------------
    # 8. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "groupby": group_col,
        "method": method,
        "min_logfc": min_logfc,
        "max_pval_adj": max_pval_adj,
        "n_genes": n_genes,
        "n_groups": len(groups),
        "results_keys": groups,
        "scanpy_version": sc.__version__,
        "omicsage_module": "deg",
        "timestamp": datetime.now().isoformat(),
    }

    adata.uns["omicsage_deg"] = provenance

    # ------------------------------------------------------------------
    # 9. Assemble deg_dict
    # ------------------------------------------------------------------
    deg_dict = {
        "results": results,
        "summary_df": summary_df,
        "provenance": provenance,
        "pairwise": pairwise_results,
    }

    return adata, deg_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _extract_group_results(
    adata: AnnData,
    group: str,
    n_genes: int,
    key: str = "rank_genes_groups",
) -> pd.DataFrame:
    """
    Extract rank_genes_groups results for a single group into a tidy DataFrame.

    Returns columns: gene, score, pval, logfc, pval_adj
    """
    rgg = adata.uns[key]

    names     = _safe_extract(rgg, "names",          group, n_genes)
    scores    = _safe_extract(rgg, "scores",         group, n_genes)
    pvals     = _safe_extract(rgg, "pvals",          group, n_genes)
    logfcs    = _safe_extract(rgg, "logfoldchanges", group, n_genes)
    pvals_adj = _safe_extract(rgg, "pvals_adj",      group, n_genes)

    df = pd.DataFrame({
        "gene":     names,
        "score":    scores,
        "pval":     pvals,
        "logfc":    logfcs,
        "pval_adj": pvals_adj,
    })

    # Drop rows with NaN gene names (scanpy pads with NaN for small groups)
    df = df.dropna(subset=["gene"])
    df = df.reset_index(drop=True)

    return df


def _safe_extract(
    rgg: dict,
    field: str,
    group: str,
    n_genes: int,
) -> np.ndarray:
    """
    Safely extract a field from rank_genes_groups structured array.
    Returns zeros/empty strings if field is missing.
    """
    if field not in rgg:
        if field in ("names",):
            return np.array([""] * n_genes)
        return np.zeros(n_genes)

    arr = rgg[field]

    # Structured array (dtype with named fields)
    if arr.dtype.names and group in arr.dtype.names:
        return arr[group]

    # Fallback — return zeros
    if field in ("names",):
        return np.array([""] * n_genes)
    return np.zeros(n_genes)


def _apply_thresholds(
    df: pd.DataFrame,
    min_logfc: float,
    max_pval_adj: float,
) -> pd.DataFrame:
    """
    Filter DEG DataFrame by log fold-change and adjusted p-value thresholds.
    Returns a copy of the filtered DataFrame, sorted by pval_adj ascending.
    """
    mask = (
        (df["pval_adj"] <= max_pval_adj) &
        (df["logfc"].abs() >= min_logfc)
    )
    filtered = df.loc[mask].copy()
    filtered = filtered.sort_values("pval_adj", ascending=True)
    filtered = filtered.reset_index(drop=True)
    return filtered


def _build_summary_df(
    results: dict[str, pd.DataFrame],
    n_top: int = 5,
) -> pd.DataFrame:
    """
    Build a summary DataFrame with top N DEGs per group.

    Returns a long-format DataFrame with columns:
        group, rank, gene, logfc, pval_adj
    """
    rows = []
    for group, df in results.items():
        top = df.head(n_top)
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            rows.append({
                "group":    group,
                "rank":     rank,
                "gene":     row["gene"],
                "logfc":    round(float(row["logfc"]), 4),
                "pval_adj": float(row["pval_adj"]),
            })
    return pd.DataFrame(rows)
