"""
pseudobulk_deg.py — Pseudobulk Differential Expression module for OmicSage
Phase 1, Step 10

Aggregates raw counts per (cell_type, donor) into bulk-like matrices and
runs pydeseq2 one-vs-rest per cell type.  Output format is identical to
deg.py deg_dict so deg_report.py works unchanged.

Usage:
    from pipeline.modules.downstream.pseudobulk_deg import pseudobulk_deg

    adata_pb, pb_dict = pseudobulk_deg(
        adata_annotated,
        groupby="cell_type_vote",
        donor_key="batch",
        min_cells=10,
        min_samples=3,
        min_logfc=0.25,
        max_pval_adj=0.05,
        inplace=False,
    )
"""

from __future__ import annotations

import warnings
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

# pydeseq2 is an optional hard-dep for this module only
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    _PYDESEQ2_AVAILABLE = True
except ImportError:
    _PYDESEQ2_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def pseudobulk_deg(
    adata: AnnData,
    groupby: str = "cell_type_vote",
    donor_key: str = "batch",
    counts_layer: str = "counts",
    min_cells: int = 10,
    min_samples: int = 3,
    min_logfc: float = 0.25,
    max_pval_adj: float = 0.05,
    n_top_genes: Optional[int] = None,
    exclude_gene_prefixes: Optional[list[str]] = None,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Run pseudobulk DEG analysis using pydeseq2, one-vs-rest per cell type.

    Aggregates raw integer counts from layers[counts_layer] per
    (groupby, donor_key) pseudo-sample, then runs DESeq2 normalisation and
    Wald-test statistics via pydeseq2.  Cell types with fewer than
    min_samples donors are skipped with a warning rather than crashing.

    Parameters
    ----------
    adata : AnnData
        Annotated AnnData.  Must contain:
          - obs[groupby]      — cell-type labels
          - obs[donor_key]    — donor / batch labels
          - layers[counts_layer] — raw integer counts (NOT log-normalised)
    groupby : str
        obs column for cell-type grouping (default: 'cell_type_vote').
    donor_key : str
        obs column for donor / sample identity (default: 'batch').
        Each unique value becomes one pseudo-bulk sample.
    counts_layer : str
        Layer holding raw integer counts (default: 'counts').
        pseudobulk sums these counts across cells in each (cell_type, donor)
        group.  Do NOT pass 'logcounts' — log-normalised values cannot be
        summed into valid bulk counts.
    min_cells : int
        Minimum number of cells required in a (cell_type, donor) combination
        to include that pseudo-sample.  Combinations below this threshold are
        dropped before DESeq2 to avoid noisy or empty pseudo-samples.
    min_samples : int
        Minimum number of donors a cell type must have (after min_cells
        filtering) to be included in the analysis.  Cell types below this
        threshold are skipped with a UserWarning — DESeq2 requires replicates.
        Default: 3 (the practical minimum for a two-group Wald test).
    min_logfc : float
        Minimum absolute log2 fold-change for filtering results.
        Applied to the LFC column returned by pydeseq2.
    max_pval_adj : float
        Maximum BH-adjusted p-value (padj) for filtering results.
    n_top_genes : int
        Maximum number of significant genes to retain per group after
        thresholding, sorted by padj ascending.  Use None to keep all.
    exclude_gene_prefixes : list of str, optional
        Gene name prefixes to remove from results after thresholding
        (e.g. ['RPL', 'RPS', 'MT-']).  Same semantics as deg.py.
    inplace : bool
        If False (default), operates on a copy of adata.

    Returns
    -------
    adata : AnnData
        AnnData with adata.uns['omicsage_pseudobulk_deg'] populated.
    pb_dict : dict
        Dictionary with keys matching deg.py deg_dict:
          'results'    — {group: DataFrame(gene, score, pval, logfc, pval_adj)}
          'summary_df' — long DataFrame: top 5 genes per group
          'provenance' — same metadata as uns
          'pairwise'   — {} (not implemented for pseudobulk; kept for compat)
          'skipped'    — {group: reason} for cell types that were skipped

    Raises
    ------
    ImportError
        If pydeseq2 is not installed.
    ValueError
        If required obs columns or counts layer are missing.
    """
    if not _PYDESEQ2_AVAILABLE:
        raise ImportError(
            "pydeseq2 is required for pseudobulk_deg().\n"
            "Install with:  pip install pydeseq2"
        )

    if not inplace:
        adata = adata.copy()

    # Normalise prefix list
    if exclude_gene_prefixes is None:
        exclude_gene_prefixes = []
    exclude_gene_prefixes = [p.upper() for p in exclude_gene_prefixes]

    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    _validate_inputs(adata, groupby, donor_key, counts_layer)

    # Cast donor and group to plain str to avoid Categorical / MultiIndex issues
    group_labels  = adata.obs[groupby].astype(str)
    donor_labels  = adata.obs[donor_key].astype(str)
    cell_types    = sorted(group_labels.unique())
    gene_names    = np.array(adata.var_names)

    # ------------------------------------------------------------------
    # 2. Extract raw counts matrix (cells × genes, dense numpy)
    # ------------------------------------------------------------------
    counts_matrix = _get_counts_matrix(adata, counts_layer)

    # ------------------------------------------------------------------
    # 3. Build pseudo-bulk count matrix: one row per (cell_type, donor)
    # ------------------------------------------------------------------
    pb_df, cell_counts_df = _aggregate_pseudobulk(
        counts_matrix, group_labels, donor_labels, gene_names, min_cells
    )

    # ------------------------------------------------------------------
    # 4. Run DESeq2 one-vs-rest per cell type
    # ------------------------------------------------------------------
    results  = {}
    skipped  = {}

    for target_ct in cell_types:
        ct_result, skip_reason = _run_deseq2_one_vs_rest(
            pb_df=pb_df,
            target_ct=target_ct,
            min_samples=min_samples,
            min_logfc=min_logfc,
            max_pval_adj=max_pval_adj,
            n_top_genes=n_top_genes,
            exclude_gene_prefixes=exclude_gene_prefixes,
        )
        if skip_reason is not None:
            skipped[target_ct] = skip_reason
            warnings.warn(
                f"pseudobulk_deg: skipping '{target_ct}' — {skip_reason}",
                UserWarning,
                stacklevel=2,
            )
        else:
            results[target_ct] = ct_result

    # ------------------------------------------------------------------
    # 5. Summary DataFrame (top 5 per group, long format)
    # ------------------------------------------------------------------
    summary_df = _build_summary_df(results, n_top=5)

    # ------------------------------------------------------------------
    # 6. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "groupby":               groupby,
        "donor_key":             donor_key,
        "counts_layer":          counts_layer,
        "method":                "DESeq2 Wald test (pydeseq2)",
        "min_cells":             min_cells,
        "min_samples":           min_samples,
        "min_logfc":             min_logfc,
        "max_pval_adj":          max_pval_adj,
        "n_top_genes":           n_top_genes,
        "exclude_gene_prefixes": exclude_gene_prefixes,
        "n_groups":              len(results),
        "n_skipped":             len(skipped),
        "results_keys":          list(results.keys()),
        "skipped_keys":          list(skipped.keys()),
        "omicsage_module":       "pseudobulk_deg",
        "timestamp":             datetime.now().isoformat(),
    }
    adata.uns["omicsage_pseudobulk_deg"] = provenance

    # ------------------------------------------------------------------
    # 7. Assemble pb_dict  (matches deg.py deg_dict layout)
    # ------------------------------------------------------------------
    pb_dict = {
        "results":    results,
        "summary_df": summary_df,
        "provenance": provenance,
        "pairwise":   {},       # not implemented; kept for deg_report compat
        "skipped":    skipped,
    }

    return adata, pb_dict


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def _validate_inputs(
    adata: AnnData,
    groupby: str,
    donor_key: str,
    counts_layer: str,
) -> None:
    """Raise ValueError with clear messages if required inputs are missing."""
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"obs column '{groupby}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )
    if donor_key not in adata.obs.columns:
        raise ValueError(
            f"obs column '{donor_key}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )
    if counts_layer not in adata.layers:
        raise ValueError(
            f"layers['{counts_layer}'] not found. "
            f"pseudobulk_deg() requires raw integer counts — "
            f"do NOT use 'logcounts'. Available layers: {list(adata.layers.keys())}"
        )


def _get_counts_matrix(adata: AnnData, counts_layer: str) -> np.ndarray:
    """Return counts as a dense float64 numpy array (cells × genes)."""
    mat = adata.layers[counts_layer]
    # Handle sparse matrices
    if hasattr(mat, "toarray"):
        mat = mat.toarray()
    return np.asarray(mat, dtype=np.float64)


def _aggregate_pseudobulk(
    counts_matrix: np.ndarray,
    group_labels: pd.Series,
    donor_labels: pd.Series,
    gene_names: np.ndarray,
    min_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Sum counts across cells in each (cell_type, donor) combination.

    Returns
    -------
    pb_df : DataFrame
        Index = (cell_type, donor) MultiIndex.
        Columns = gene names.
        Values = summed raw counts (float64).
    cell_counts_df : DataFrame
        Index = (cell_type, donor) MultiIndex.
        Column = 'n_cells' — number of cells that contributed to each pseudo-sample.
        Combinations with n_cells < min_cells are already removed.
    """
    # Build a combined grouping series
    combo = group_labels.str.cat(donor_labels, sep="||")
    unique_combos = sorted(combo.unique())

    rows    = []
    n_cells = []
    index   = []

    for c in unique_combos:
        mask = (combo == c).values
        n    = int(mask.sum())
        if n < min_cells:
            continue
        cell_type, donor = c.split("||", 1)
        summed = counts_matrix[mask].sum(axis=0)   # shape (n_genes,)
        rows.append(summed)
        n_cells.append(n)
        index.append((cell_type, donor))

    if not rows:
        raise ValueError(
            f"No (cell_type, donor) combinations survived the min_cells={min_cells} "
            f"filter.  Check that obs columns exist and min_cells is not too large."
        )

    mi = pd.MultiIndex.from_tuples(index, names=["cell_type", "donor"])
    pb_df           = pd.DataFrame(np.vstack(rows), index=mi, columns=gene_names)
    cell_counts_df  = pd.DataFrame({"n_cells": n_cells}, index=mi)

    return pb_df, cell_counts_df


# ---------------------------------------------------------------------------
# DESeq2 one-vs-rest runner
# ---------------------------------------------------------------------------

def _run_deseq2_one_vs_rest(
    pb_df: pd.DataFrame,
    target_ct: str,
    min_samples: int,
    min_logfc: float,
    max_pval_adj: float,
    n_top_genes: Optional[int],
    exclude_gene_prefixes: list[str],
) -> tuple[Optional[pd.DataFrame], Optional[str]]:
    """
    Run a single DESeq2 Wald test: target_ct vs all other cell types.

    Returns (result_df, None) on success or (None, reason_string) if skipped.
    result_df columns: gene, score (stat), pval, logfc (log2FoldChange), pval_adj
    """
    # Subset pseudo-bulk rows for target vs rest
    is_target = pb_df.index.get_level_values("cell_type") == target_ct
    n_target  = int(is_target.sum())
    n_rest    = int((~is_target).sum())

    if n_target < min_samples:
        return None, (
            f"{n_target} donor pseudo-samples (need ≥ {min_samples}). "
            f"Increase donors or lower min_samples."
        )
    if n_rest < min_samples:
        return None, (
            f"rest group has only {n_rest} donor pseudo-samples (need ≥ {min_samples})."
        )

    # Build count matrix and metadata for pydeseq2
    # pydeseq2 expects: counts (samples × genes, int), metadata (samples × covariates)
    subset_df  = pb_df.copy()
    condition  = np.where(is_target, target_ct, "rest")

    counts_int = subset_df.values.astype(int)   # DESeq2 requires integers

    # Remove genes with zero counts across all samples in this comparison
    gene_mask   = counts_int.sum(axis=0) > 0
    counts_int  = counts_int[:, gene_mask]
    gene_names  = np.array(subset_df.columns)[gene_mask]

    sample_names = [f"s{i}" for i in range(len(condition))]

    counts_df  = pd.DataFrame(
        counts_int,
        index=sample_names,
        columns=gene_names,
    )
    metadata_df = pd.DataFrame(
        {"condition": condition},
        index=sample_names,
    )

    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design_factors="condition",
            quiet=True,
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=["condition", target_ct, "rest"], quiet=True)
        stat_res.summary()
        try:
            stat_res.lfc_shrink(coeff=f"condition_{target_ct}_vs_rest")
        except Exception:
            pass  # shrinkage may fail on some edge cases — fall back to MLE LFC

        res_df = stat_res.results_df.copy()
    except Exception as e:
        return None, f"pydeseq2 error: {e}"

    # Rename columns to match deg.py output schema
    # pydeseq2 returns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    rename_map = {
        "log2FoldChange": "logfc",
        "pvalue":         "pval",
        "padj":           "pval_adj",
        "stat":           "score",
    }
    res_df = res_df.rename(columns=rename_map)
    res_df = res_df.reset_index().rename(columns={"index": "gene"})

    # Keep only the columns deg.py uses
    keep_cols = ["gene", "score", "pval", "logfc", "pval_adj"]
    for col in keep_cols:
        if col not in res_df.columns:
            res_df[col] = np.nan
    res_df = res_df[keep_cols].copy()

    # Drop genes where pydeseq2 returned NaN padj (low-count outliers)
    res_df = res_df.dropna(subset=["pval_adj"]).copy()

    # Apply thresholds
    res_df = _apply_thresholds(res_df, min_logfc, max_pval_adj)

    # Apply gene prefix exclusion
    res_df = _apply_prefix_exclusion(res_df, exclude_gene_prefixes, target_ct)

    # Cap at n_top_genes
    if n_top_genes is not None:
        res_df = res_df.head(n_top_genes)

    res_df = res_df.reset_index(drop=True)
    return res_df, None


# ---------------------------------------------------------------------------
# Shared helpers (mirrors deg.py style)
# ---------------------------------------------------------------------------

def _apply_thresholds(
    df: pd.DataFrame,
    min_logfc: float,
    max_pval_adj: float,
) -> pd.DataFrame:
    """Filter by |logfc| and pval_adj; sort by pval_adj ascending."""
    mask = (
        (df["pval_adj"] <= max_pval_adj) &
        (df["logfc"].abs() >= min_logfc)
    )
    filtered = df.loc[mask].copy()
    filtered = filtered.sort_values("pval_adj", ascending=True)
    return filtered.reset_index(drop=True)


def _apply_prefix_exclusion(
    df: pd.DataFrame,
    exclude_prefixes: list[str],
    group: str,
) -> pd.DataFrame:
    """Remove genes whose names start with any of the given prefixes."""
    if not exclude_prefixes or df.empty:
        return df

    mask_keep = ~df["gene"].str.upper().apply(
        lambda g: any(g.startswith(p) for p in exclude_prefixes)
    )
    n_removed = int((~mask_keep).sum())

    if n_removed > 0:
        removed = df.loc[~mask_keep, "gene"].tolist()
        warnings.warn(
            f"Group '{group}': excluded {n_removed} gene(s) matching prefixes "
            f"{exclude_prefixes}: {removed[:10]}"
            f"{'...' if len(removed) > 10 else ''}",
            UserWarning,
            stacklevel=3,
        )

    return df.loc[mask_keep].copy().reset_index(drop=True)


def _build_summary_df(
    results: dict[str, pd.DataFrame],
    n_top: int = 5,
) -> pd.DataFrame:
    """
    Build a long-format summary DataFrame with top N DEGs per group.
    Columns: group, rank, gene, logfc, pval_adj
    Matches deg.py _build_summary_df output exactly so deg_report.py works.
    """
    rows = []
    for group, df in results.items():
        for rank, (_, row) in enumerate(df.head(n_top).iterrows(), start=1):
            rows.append({
                "group":    group,
                "rank":     rank,
                "gene":     row["gene"],
                "logfc":    round(float(row["logfc"]), 4),
                "pval_adj": float(row["pval_adj"]),
            })
    return pd.DataFrame(rows)
