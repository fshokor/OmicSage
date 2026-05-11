"""
gsea.py — Gene Set Enrichment Analysis module for OmicSage
Phase 1, Step 7

Performs over-representation analysis (ORA) via gseapy.enrichr.
Accepts deg_dict['results'] from deg() and runs per-group enrichment
against GO Biological Process, KEGG, and Reactome (configurable).

Usage:
    from pipeline.modules.downstream.gsea import gsea

    adata_gsea, gsea_dict = gsea(
        adata_deg,
        deg_dict=deg_dict,
        gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"],
        min_logfc=0.25,
        max_pval_adj=0.05,
        top_n_genes=None,
        min_genes=5,
        organism="human",
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

try:
    import gseapy
    _GSEAPY_AVAILABLE = True
except ImportError:
    _GSEAPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

_DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
]

_REQUIRED_ENRICHR_COLS = {
    "Term", "Overlap", "P-value", "Adjusted P-value", "Genes",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def gsea(
    adata: AnnData,
    deg_dict: dict,
    gene_sets: Optional[list[str]] = None,
    min_logfc: float = 0.25,
    max_pval_adj: float = 0.05,
    top_n_genes: Optional[int] = None,
    min_genes: int = 5,
    organism: str = "human",
    direction: str = "up",
    exclude_gene_prefixes: Optional[list[str]] = None,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Run over-representation analysis (ORA) on DEG results using gseapy.enrichr.

    Parameters
    ----------
    adata : AnnData
        AnnData returned by deg(). Used to derive the gene universe
        (all detected genes in adata.var_names).
    deg_dict : dict
        deg_dict returned by deg(). Must contain 'results' key:
        {group: DataFrame(gene, score, pval, logfc, pval_adj)}.
    gene_sets : list of str, optional
        Enrichr gene set library names to query.
        Defaults to GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022.
        Any valid Enrichr library name is accepted — e.g. "DisGeNET",
        "KEGG_2021_Mouse", "ENCODE_TF_ChIP-seq_2015".
        Run gseapy.get_library_name() to list all available names.
    min_logfc : float
        Minimum absolute log2 fold-change for a DEG to be included.
        Applied in addition to any filtering already done by deg().
    max_pval_adj : float
        Maximum adjusted p-value for a DEG to be included in the query list.
    top_n_genes : int or None
        If set, use only the top N DEGs per group (ranked by pval_adj).
        If None, use all DEGs passing filters.
    min_genes : int
        Minimum number of query genes required to run enrichment for a group.
        Groups with fewer genes are skipped with a UserWarning.
    organism : str
        Organism for Enrichr query. One of: "human", "mouse", "yeast", etc.
        Note: gene set names must match the organism.
    direction : str
        Which DEGs to use as the ORA query gene list. One of:

        "up"   — upregulated genes only (logfc >= +min_logfc).
                 Default. Best for cell-type marker enrichment.

        "down" — downregulated genes only (logfc <= -min_logfc).
                 Use for suppressed pathways, e.g. tumour suppressor
                 loss in cancer data.

        "both" — runs two separate ORA queries per group: one for
                 upregulated and one for downregulated genes.
                 Results are stored under keys "{group}__up" and
                 "{group}__down" in gsea_dict['results'].
                 The summary_df gains a 'Direction' column.
                 Use when you want a complete pathway picture, e.g.
                 comparing cancer vs normal where both activation and
                 suppression of pathways are biologically meaningful.

        Note: up and down queries are always run independently — mixing
        both directions into a single ORA query would produce uninterpretable
        results because a pathway could match genes going in opposite
        directions simultaneously.
    exclude_gene_prefixes : list of str, optional
        Gene name prefixes to remove from the query gene list before running
        enrichment. Applied on top of any filtering already done by deg().
        Example: ["RPL", "RPS", "MT-"]
        Note: these genes remain in the gene universe (background) — only
        the query list is filtered, which is statistically correct for ORA.
        Default: None (no additional filtering).
    inplace : bool
        If False (default), operates on a copy of adata.
        If True, modifies adata in place.

    Returns
    -------
    adata : AnnData
        AnnData with adata.uns['omicsage_gsea'] populated.
    gsea_dict : dict
        Dictionary with keys:
            'results'    — {group: DataFrame} for direction="up"/"down", or
                           {group__up: DataFrame, group__down: DataFrame}
                           for direction="both"
            'summary_df' — top 3 terms per group (+ Direction column if "both")
            'provenance' — same metadata as uns['omicsage_gsea']
            'skipped'    — list of (group, direction) tuples skipped due to
                           fewer than min_genes DEGs

    Raises
    ------
    ImportError
        If gseapy is not installed.
    ValueError
        If deg_dict does not contain a 'results' key, or direction is invalid.
    """
    if not _GSEAPY_AVAILABLE:
        raise ImportError(
            "gseapy is required for GSEA. "
            "Install with: pip install gseapy  "
            "or: conda install -c bioconda gseapy"
        )

    if "results" not in deg_dict:
        raise ValueError(
            "deg_dict must contain a 'results' key. "
            "Pass the deg_dict returned by deg()."
        )

    valid_directions = {"up", "down", "both"}
    if direction not in valid_directions:
        raise ValueError(
            f"direction must be one of {valid_directions}, got '{direction}'."
        )

    if not inplace:
        adata = adata.copy()

    if gene_sets is None:
        gene_sets = _DEFAULT_GENE_SETS

    # Ensure gene_sets is always a list
    if isinstance(gene_sets, str):
        gene_sets = [gene_sets]

    # Normalise prefix list
    if exclude_gene_prefixes is None:
        exclude_gene_prefixes = []
    exclude_gene_prefixes = [p.upper() for p in exclude_gene_prefixes]

    # ------------------------------------------------------------------
    # 1. Validate requested gene sets against Enrichr
    # ------------------------------------------------------------------
    _validate_gene_sets(gene_sets)

    # ------------------------------------------------------------------
    # 2. Build gene universe from adata.var_names
    # ------------------------------------------------------------------
    gene_universe = adata.var_names.tolist()

    # ------------------------------------------------------------------
    # 3. Determine which directions to run
    # ------------------------------------------------------------------
    directions_to_run = ["up", "down"] if direction == "both" else [direction]

    # ------------------------------------------------------------------
    # 4. Run ORA per group × direction
    # ------------------------------------------------------------------
    deg_results = deg_dict["results"]
    results: dict[str, pd.DataFrame] = {}
    skipped: list[tuple[str, str]] = []

    for group, df in deg_results.items():
        for dir_ in directions_to_run:
            query_genes = _build_query_genes(
                df,
                min_logfc=min_logfc,
                max_pval_adj=max_pval_adj,
                top_n=top_n_genes,
                direction=dir_,
                exclude_prefixes=exclude_gene_prefixes,
            )

            # Key is plain group name for single direction,
            # suffixed for "both" so up and down are stored separately
            result_key = f"{group}__{dir_}" if direction == "both" else group

            if len(query_genes) < min_genes:
                warnings.warn(
                    f"Group '{group}' ({dir_}): only {len(query_genes)} DEG(s) "
                    f"passing filters (min_genes={min_genes}). Skipping.",
                    UserWarning,
                    stacklevel=2,
                )
                skipped.append((group, dir_))
                continue

            group_df = _run_enrichr(
                query_genes=query_genes,
                gene_sets=gene_sets,
                gene_universe=gene_universe,
                organism=organism,
                group=result_key,
            )

            # Tag each row with direction so the report can colour-code
            group_df["Direction"] = dir_
            results[result_key] = group_df

    # ------------------------------------------------------------------
    # 5. Build summary DataFrame (top 3 terms per group × direction)
    # ------------------------------------------------------------------
    summary_df = _build_summary_df(results, n_top=3, direction_mode=direction)

    # ------------------------------------------------------------------
    # 6. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "groupby": deg_dict.get("provenance", {}).get("groupby", "unknown"),
        "gene_sets": gene_sets,
        "organism": organism,
        "direction": direction,
        "min_logfc": min_logfc,
        "max_pval_adj": max_pval_adj,
        "top_n_genes": top_n_genes,
        "min_genes": min_genes,
        "exclude_gene_prefixes": exclude_gene_prefixes,
        "n_groups_tested": len(results),
        "n_groups_skipped": len(skipped),
        "groups_skipped": skipped,
        "gseapy_version": gseapy.__version__ if _GSEAPY_AVAILABLE else "unavailable",
        "omicsage_module": "gsea",
        "timestamp": datetime.now().isoformat(),
    }

    adata.uns["omicsage_gsea"] = provenance

    # ------------------------------------------------------------------
    # 7. Assemble gsea_dict
    # ------------------------------------------------------------------
    gsea_dict = {
        "results": results,
        "summary_df": summary_df,
        "provenance": provenance,
        "skipped": skipped,
    }

    return adata, gsea_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _validate_gene_sets(gene_sets: list[str]) -> None:
    """
    Check requested gene set names against Enrichr.
    Warns (does not raise) if a name is not found — Enrichr updates names
    periodically, so a warning is more useful than a hard crash.
    """
    try:
        available = gseapy.get_library_name(organism="human")
        available_set = set(available)
        for gs in gene_sets:
            if gs not in available_set:
                warnings.warn(
                    f"Gene set '{gs}' not found in Enrichr library list. "
                    f"It may be misspelled or have been renamed. "
                    f"Run gseapy.get_library_name() to see all available names.",
                    UserWarning,
                    stacklevel=3,
                )
    except Exception:
        # Network unavailable or Enrichr down — warn and continue
        warnings.warn(
            "Could not validate gene set names against Enrichr "
            "(network unavailable?). Proceeding anyway.",
            UserWarning,
            stacklevel=3,
        )


def _build_query_genes(
    df: pd.DataFrame,
    min_logfc: float,
    max_pval_adj: float,
    top_n: Optional[int],
    direction: str = "up",
    exclude_prefixes: Optional[list[str]] = None,
) -> list[str]:
    """
    Filter a DEG DataFrame and return a directional query gene list.

    Parameters
    ----------
    direction : "up" | "down"
        "up"   → logfc >= +min_logfc  (upregulated markers)
        "down" → logfc <= -min_logfc  (downregulated / suppressed genes)

    exclude_prefixes filters the query list only. The gene universe
    (background) is NOT affected — this is statistically correct because
    removing genes from the background would artificially inflate
    enrichment scores for pathways containing those genes.
    """
    if df.empty:
        return []

    pval_mask = df["pval_adj"] <= max_pval_adj

    if direction == "up":
        dir_mask = df["logfc"] >= min_logfc
    else:  # "down"
        dir_mask = df["logfc"] <= -min_logfc

    filtered = df.loc[pval_mask & dir_mask].copy()
    filtered = filtered.sort_values("pval_adj", ascending=True)

    # Apply prefix exclusion to query list only (not to universe)
    if exclude_prefixes:
        keep = ~filtered["gene"].str.upper().apply(
            lambda g: any(g.startswith(p) for p in exclude_prefixes)
        )
        filtered = filtered.loc[keep]

    if top_n is not None:
        filtered = filtered.head(top_n)

    return filtered["gene"].dropna().tolist()


def _run_enrichr(
    query_genes: list[str],
    gene_sets: list[str],
    gene_universe: list[str],
    organism: str,
    group: str,
) -> pd.DataFrame:
    """
    Run gseapy.enrichr for one group and return a tidy DataFrame.

    Returns columns: Term, Overlap, P-value, Adjusted P-value, Genes
    Falls back to an empty DataFrame if enrichr fails.
    """
    try:
        enr = gseapy.enrichr(
            gene_list=query_genes,
            gene_sets=gene_sets,
            background=gene_universe,
            organism=organism,
            outdir=None,        # do not write files to disk
            verbose=False,
        )
        df = enr.results.copy()

        # Standardise column names — gseapy uses slightly different names
        # across versions; normalise to our expected set
        df = _normalise_enrichr_columns(df)

        # Sort by adjusted p-value
        if "Adjusted P-value" in df.columns:
            df = df.sort_values("Adjusted P-value", ascending=True)
            df = df.reset_index(drop=True)

        return df

    except Exception as e:
        warnings.warn(
            f"enrichr failed for group '{group}': {e}. "
            f"Returning empty DataFrame for this group.",
            UserWarning,
            stacklevel=3,
        )
        return pd.DataFrame(columns=list(_REQUIRED_ENRICHR_COLS))


def _normalise_enrichr_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise gseapy.enrichr output column names across versions.

    gseapy ≥1.0 (confirmed on 1.1.3) returns these columns:
      Gene_set, Term, P-value, Adjusted P-value, Old P-value,
      Old adjusted P-value, Odds Ratio, Combined Score, Genes

    Notably: NO "Overlap" column in gseapy ≥1.0.
    We derive it by counting the semicolon-separated genes in "Genes".

    Older versions may use:
      "Adj. P-value"     → "Adjusted P-value"
      "Genes_in_Overlap" → "Genes"
    """
    rename_map = {
        "Adj. P-value":     "Adjusted P-value",
        "Genes_in_Overlap": "Genes",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    # Derive Overlap from Genes column if absent or all-NaN
    # gseapy ≥1.0 dropped the Overlap column entirely.
    # We reconstruct it as "k matched genes" (denominator unknown without
    # the full gene set size, so we just report the count).
    if "Overlap" not in df.columns or df["Overlap"].isna().all():
        if "Genes" in df.columns:
            def _count_genes(genes_str) -> str:
                if not isinstance(genes_str, str) or not genes_str.strip():
                    return "0"
                return str(len([g for g in genes_str.split(";") if g.strip()]))
            df["Overlap"] = df["Genes"].apply(_count_genes)
        else:
            df["Overlap"] = "0"

    # Ensure all required columns are present
    for col in _REQUIRED_ENRICHR_COLS:
        if col not in df.columns:
            df[col] = np.nan

    return df


def _build_summary_df(
    results: dict[str, pd.DataFrame],
    n_top: int = 3,
    direction_mode: str = "up",
) -> pd.DataFrame:
    """
    Build a long-format summary DataFrame with top N terms per result key.

    When direction_mode is "both", result keys are "{group}__up" /
    "{group}__down" and the summary gains a 'Direction' column and a
    clean 'group' column (without the suffix) for display.

    Returns columns: group, direction (if both), rank, Term, Overlap,
                     Adjusted P-value, Genes
    """
    rows = []
    for key, df in results.items():
        if df.empty:
            continue

        # Parse group name and direction from the result key
        if direction_mode == "both" and "__" in key:
            group, dir_label = key.rsplit("__", 1)
        else:
            group = key
            dir_label = direction_mode

        top = df.head(n_top)
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            entry = {
                "group":            group,
                "direction":        dir_label,
                "rank":             rank,
                "Term":             row.get("Term", ""),
                "Overlap":          row.get("Overlap", ""),
                "Adjusted P-value": row.get("Adjusted P-value", np.nan),
                "Genes":            row.get("Genes", ""),
            }
            rows.append(entry)

    return pd.DataFrame(rows)
