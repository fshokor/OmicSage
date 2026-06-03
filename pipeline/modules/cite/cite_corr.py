"""
cite_corr.py — Protein-RNA Correlation module for OmicSage
Phase 4, Step 9

Computes per-protein Spearman correlation between ADT expression and the
corresponding gene's RNA expression, computed WITHIN each cell type and
then aggregated using Fisher z-transform averaging.

Why within-cell-type?
---------------------
Computing Spearman r across all cells mixes cell populations. This creates
systematic artifacts when the panel is dominated by markers of one lineage
(e.g. T cell markers in BMMC data): CLR normalization suppresses all other
proteins in T cells relative to the population mean, producing spurious
negative correlations. Within-cell-type computation isolates the
post-transcriptional regulation signal from the compositional contrast
between cell types.

Fisher z-transform aggregation
-------------------------------
Per-cell-type r values are combined using the Fisher z-transform:
  z_i = arctanh(r_i)
  z_mean = weighted mean of z_i  (weight = sqrt(n_i - 3))
  r_agg = tanh(z_mean)
This is the statistically correct way to average correlations across groups
of different sizes.

Biological interpretation
-------------------------
High aggregate r (> 0.4):
    Protein and mRNA co-regulated within cell types — typical for
    housekeeping or lineage markers measured at steady state.

Low or zero aggregate r (|r| < 0.2):
    Post-transcriptional regulation — protein abundance controlled at
    translation, secretion, or membrane trafficking. Biologically interesting.
    Expected for CD4, CD8, CD3 complex subunits.

Negative aggregate r:
    True within-cell-type anti-correlation. Rare and biologically meaningful —
    may indicate feedback inhibition or splice variants. Distinguishable from
    the cross-cell-type compositional artifact fixed by this module.

Input
-----
MuData from cite_06 integration step (preferred).
Requires:
  mdata["adt"].layers["adt_clr"]   — CLR-normalised ADT (float matrix)
  mdata["rna"].X or layers["logcounts"]  — log-normalised RNA
  mdata["adt"].var_names            — protein names (e.g. "CD3E", "CD19")
  mdata["rna"].var_names            — gene symbols (e.g. "CD3E", "CD19")
  cell type label in obs            — resolved via celltype_key parameter

Outputs written to data.uns
---------------------------
uns["omicsage_cite_corr"]            Provenance block
uns["omicsage_cite_corr_results"]    Aggregate results as JSON string
uns["omicsage_cite_corr_per_ct"]     Per-cell-type results as JSON string

Returns
-------
(MuData, corr_dict)

corr_dict keys
--------------
  "results"              : DataFrame(protein, gene, r, pval, pval_adj,
                                     n_cells, n_celltypes, matched_by)
                           r is the Fisher-z-aggregated value across cell types
  "results_per_celltype" : DataFrame(protein, gene, cell_type, r, n_cells)
  "matched"              : list[tuple[str, str]] — (protein, gene) pairs tested
  "unmatched"            : list[str] — protein names with no RNA match
  "skipped"              : list[dict] — pairs skipped (below min_cells globally)
  "provenance"           : same as uns["omicsage_cite_corr"]

Usage
-----
    from pipeline.modules.cite.cite_corr import cite_corr

    mdata_corr, corr_dict = cite_corr(
        mdata,
        method="spearman",
        min_cells=30,
        min_cells_per_ct=10,
        celltype_key="auto",
        use_logcounts=True,
        inplace=False,
    )
"""

from __future__ import annotations

import re
import warnings
from datetime import datetime, timezone
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from anndata import AnnData

try:
    from mudata import MuData
    _MUDATA_AVAILABLE = True
except ImportError:
    _MUDATA_AVAILABLE = False
    MuData = None  # type: ignore

try:
    from statsmodels.stats.multitest import multipletests
    _STATSMODELS_AVAILABLE = True
except ImportError:
    _STATSMODELS_AVAILABLE = False


# Cell type obs columns to try, in priority order
_CELLTYPE_KEY_CANDIDATES = [
    "adt_celltype_manual",
    "adt_celltype",
    "cell_type_vote",
    "cell_type",
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cite_corr(
    data: Union["MuData", AnnData],
    method: str = "spearman",
    min_cells: int = 30,
    min_cells_per_ct: int = 10,
    celltype_key: str = "auto",
    use_logcounts: bool = True,
    inplace: bool = False,
) -> tuple[Union["MuData", AnnData], dict]:
    """
    Compute per-protein Spearman correlation between ADT and RNA expression,
    within each cell type, aggregated via Fisher z-transform.

    Parameters
    ----------
    data : MuData
        MuData from cite_06 integration step.
        Requires mdata["adt"] and mdata["rna"] modalities.
        AnnData input raises ValueError immediately.
    method : str
        Correlation method. Only "spearman" is currently supported.
    min_cells : int
        Minimum number of cells expressing both the protein (ADT CLR > 0)
        and the gene (RNA > 0) across the *entire dataset* for a pair to be
        considered at all. Pairs below this threshold are skipped entirely.
        Default: 30.
    min_cells_per_ct : int
        Minimum number of cells expressing both modalities *within a single
        cell type* for that cell type to contribute to the Fisher-z aggregate.
        Cell types below this threshold are excluded from aggregation for that
        pair (but other cell types still contribute).
        Default: 10.
    celltype_key : str
        obs column to use for cell type grouping. "auto" (default) tries
        the following in order:
          1. mdata["adt"].obs["adt_celltype_manual"]
          2. mdata["adt"].obs["adt_celltype"]
          3. mdata["rna"].obs["cell_type_vote"]
          4. mdata["rna"].obs["cell_type"]
        If none found, falls back to global correlation with a UserWarning.
    use_logcounts : bool
        If True (default), uses mdata["rna"].layers["logcounts"] for RNA.
        Falls back to mdata["rna"].X with a UserWarning if absent.
    inplace : bool
        If False (default), operates on a copy.

    Returns
    -------
    data : MuData
        Input MuData with uns["omicsage_cite_corr"] populated.
    corr_dict : dict
        Keys: results, results_per_celltype, matched, unmatched,
              skipped, provenance.
    """
    # ------------------------------------------------------------------
    # 0. Validate input type
    # ------------------------------------------------------------------
    if not (_MUDATA_AVAILABLE and isinstance(data, MuData)):
        raise ValueError(
            "cite_corr requires a MuData object with both 'rna' and 'adt' "
            "modalities (from cite_06 integration step). "
            "AnnData (single-modality) input is not supported."
        )

    if method != "spearman":
        raise ValueError(
            f"method='{method}' is not supported. "
            "Only 'spearman' is currently implemented."
        )

    if not inplace:
        data = data.copy()

    rna = data["rna"]
    adt = data["adt"]

    # ------------------------------------------------------------------
    # 1. Resolve cell type labels
    # ------------------------------------------------------------------
    cell_types, resolved_key = _resolve_celltype_key(adt, rna, celltype_key)
    use_per_celltype = cell_types is not None

    if not use_per_celltype:
        warnings.warn(
            "cite_corr: no cell type label column found. "
            "Falling back to global correlation (may produce artifacts). "
            f"Tried: {_CELLTYPE_KEY_CANDIDATES}. "
            "Set celltype_key explicitly or run adt_annotate first.",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 2. Extract expression matrices
    # ------------------------------------------------------------------
    adt_matrix = _get_adt_matrix(adt)          # (n_cells, n_proteins)
    rna_matrix = _get_rna_matrix(rna, use_logcounts)  # (n_cells, n_genes)

    adt_proteins = list(adt.var_names)
    rna_genes    = list(rna.var_names)
    rna_gene_lower = {g.lower(): i for i, g in enumerate(rna_genes)}

    # ------------------------------------------------------------------
    # 3. Match proteins to genes
    # ------------------------------------------------------------------
    matched:   list[tuple[str, str, int, int, str]] = []
    unmatched: list[str] = []

    for prot_idx, protein in enumerate(adt_proteins):
        gene_idx, matched_by = _match_protein_to_gene(
            protein, rna_gene_lower, rna_genes
        )
        if gene_idx is not None:
            matched.append((protein, rna_genes[gene_idx], prot_idx, gene_idx, matched_by))
        else:
            unmatched.append(protein)

    if not matched:
        warnings.warn(
            "cite_corr: no protein names matched any RNA gene symbols. "
            "Check that mdata['adt'].var_names use standard gene-symbol names.",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 4. Compute correlations (per cell type → Fisher-z aggregate)
    # ------------------------------------------------------------------
    rows: list[dict] = []           # aggregate results (one row per pair)
    per_ct_rows: list[dict] = []    # per-cell-type results
    skipped: list[dict] = []

    unique_celltypes = (
        sorted(cell_types.unique().tolist()) if use_per_celltype else [None]
    )

    for protein, gene, prot_idx, gene_idx, matched_by in matched:
        x_prot = _get_column(adt_matrix, prot_idx)
        x_rna  = _get_column(rna_matrix, gene_idx)

        # Global expression mask — need enough cells across full dataset
        global_mask = (x_prot > 0) & (x_rna > 0)
        n_cells_global = int(global_mask.sum())

        if n_cells_global < min_cells:
            skipped.append({
                "protein": protein,
                "gene":    gene,
                "n_cells": n_cells_global,
                "reason":  f"< min_cells ({min_cells}) globally",
            })
            continue

        if use_per_celltype:
            # --- Within-cell-type correlations ---
            z_values:  list[float] = []
            weights:   list[float] = []

            for ct in unique_celltypes:
                ct_mask = (cell_types == ct).values
                mask = ct_mask & (x_prot > 0) & (x_rna > 0)
                n_ct = int(mask.sum())

                if n_ct < min_cells_per_ct:
                    continue

                r_ct, _ = stats.spearmanr(x_prot[mask], x_rna[mask])

                # Clamp r to avoid arctanh blow-up at ±1
                r_ct = float(np.clip(r_ct, -0.9999, 0.9999))

                per_ct_rows.append({
                    "protein":   protein,
                    "gene":      gene,
                    "cell_type": ct,
                    "r":         r_ct,
                    "n_cells":   n_ct,
                })

                # Fisher z-transform weight = sqrt(n - 3)
                w = max(float(np.sqrt(n_ct - 3)), 0.0)
                z_values.append(np.arctanh(r_ct))
                weights.append(w)

            if not z_values:
                skipped.append({
                    "protein": protein,
                    "gene":    gene,
                    "n_cells": n_cells_global,
                    "reason":  f"no cell type had >= min_cells_per_ct ({min_cells_per_ct})",
                })
                continue

            # Weighted Fisher-z mean → aggregate r
            weights_arr = np.array(weights)
            z_arr       = np.array(z_values)
            total_w     = weights_arr.sum()
            if total_w == 0:
                r_agg = 0.0
            else:
                z_mean = np.dot(weights_arr, z_arr) / total_w
                r_agg  = float(np.tanh(z_mean))

            n_celltypes_used = len(z_values)

            # Approximate p-value: run Spearman on all cells that expressed
            # both modalities across any contributing cell type.
            # This is a conservative global p-value — it uses more cells than
            # any single cell type, giving a lower bound on significance.
            # The aggregate r (Fisher-z) is the primary statistic; p is
            # reported for completeness.
            _, pval = stats.spearmanr(x_prot[global_mask], x_rna[global_mask])

        else:
            # --- Fallback: global correlation ---
            mask = global_mask
            r_agg, pval = stats.spearmanr(x_prot[mask], x_rna[mask])
            r_agg = float(r_agg)
            n_celltypes_used = 0

        rows.append({
            "protein":       protein,
            "gene":          gene,
            "r":             r_agg,
            "pval":          float(pval),
            "n_cells":       n_cells_global,
            "n_celltypes":   n_celltypes_used,
            "matched_by":    matched_by,
        })

    # ------------------------------------------------------------------
    # 5. Multiple testing correction (BH FDR) on aggregate results
    # ------------------------------------------------------------------
    if rows:
        results_df = pd.DataFrame(rows)
        if _STATSMODELS_AVAILABLE and len(results_df) > 1:
            _, pval_adj, _, _ = multipletests(
                results_df["pval"].values, method="fdr_bh"
            )
            results_df["pval_adj"] = pval_adj
        else:
            results_df["pval_adj"] = results_df["pval"]
            if not _STATSMODELS_AVAILABLE:
                warnings.warn(
                    "statsmodels not available — BH FDR correction skipped.",
                    UserWarning, stacklevel=2,
                )
        results_df = results_df.sort_values("r", ascending=False).reset_index(drop=True)
    else:
        results_df = pd.DataFrame(
            columns=["protein", "gene", "r", "pval", "pval_adj",
                     "n_cells", "n_celltypes", "matched_by"]
        )

    per_ct_df = (
        pd.DataFrame(per_ct_rows)
        if per_ct_rows
        else pd.DataFrame(columns=["protein", "gene", "cell_type", "r", "n_cells"])
    )

    # ------------------------------------------------------------------
    # 6. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "module":              "cite_corr",
        "timestamp":           datetime.now(timezone.utc).isoformat(),
        "method":              method,
        "min_cells":           min_cells,
        "min_cells_per_ct":    min_cells_per_ct,
        "celltype_key":        resolved_key or "none (global fallback)",
        "use_per_celltype":    use_per_celltype,
        "use_logcounts":       use_logcounts,
        "n_proteins_panel":    len(adt_proteins),
        "n_matched":           len(matched),
        "n_unmatched":         len(unmatched),
        "n_tested":            len(rows),
        "n_skipped":           len(skipped),
        "n_celltypes":         len(unique_celltypes) if use_per_celltype else 0,
        "bh_correction":       _STATSMODELS_AVAILABLE and len(rows) > 1,
        "unmatched_proteins":  unmatched,
    }
    data.uns["omicsage_cite_corr"] = provenance

    data.uns["omicsage_cite_corr_results"] = (
        results_df.to_json(orient="records") if not results_df.empty else "[]"
    )
    data.uns["omicsage_cite_corr_per_ct"] = (
        per_ct_df.to_json(orient="records") if not per_ct_df.empty else "[]"
    )

    # ------------------------------------------------------------------
    # 7. Assemble return dict
    # ------------------------------------------------------------------
    corr_dict = {
        "results":              results_df,
        "results_per_celltype": per_ct_df,
        "matched":              [(p, g) for p, g, _, _, _ in matched],
        "unmatched":            unmatched,
        "skipped":              skipped,
        "provenance":           provenance,
    }

    return data, corr_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

_STRIP_SUFFIXES = re.compile(
    r"[-_](TotalSeqA|TotalSeqB|TotalSeqC|TotalSeq|"
    r"BV421|BV711|APC|PE|FITC|BioLegend|[0-9]+).*$",
    re.IGNORECASE,
)


def _resolve_celltype_key(
    adt: AnnData,
    rna: AnnData,
    celltype_key: str,
) -> tuple[Optional[pd.Series], Optional[str]]:
    """
    Resolve the cell type label Series and the key name used.

    Returns (cell_type_series, key_name) or (None, None) if not found.
    """
    if celltype_key != "auto":
        # Explicit key — check adt.obs first, then rna.obs
        for obj in (adt, rna):
            if celltype_key in obj.obs.columns:
                return obj.obs[celltype_key].astype(str), celltype_key
        warnings.warn(
            f"cite_corr: celltype_key='{celltype_key}' not found in "
            "adt.obs or rna.obs. Falling back to auto-detection.",
            UserWarning, stacklevel=3,
        )

    # Auto-detect
    for key in _CELLTYPE_KEY_CANDIDATES:
        if key in adt.obs.columns:
            return adt.obs[key].astype(str), key
        if key in rna.obs.columns:
            return rna.obs[key].astype(str), key

    return None, None


def _match_protein_to_gene(
    protein: str,
    rna_gene_lower: dict[str, int],
    rna_genes: list[str],
) -> tuple[Optional[int], Optional[str]]:
    """Try to match a protein name to a gene symbol."""
    key = protein.lower()
    if key in rna_gene_lower:
        return rna_gene_lower[key], "exact"

    stripped = _STRIP_SUFFIXES.sub("", protein).strip()
    if stripped.lower() != key and stripped.lower() in rna_gene_lower:
        return rna_gene_lower[stripped.lower()], "stripped"

    return None, None


def _get_adt_matrix(adt: AnnData) -> np.ndarray:
    """Return dense (n_cells, n_proteins) CLR matrix."""
    if "adt_clr" in adt.layers:
        mat = adt.layers["adt_clr"]
    else:
        warnings.warn(
            "adt.layers['adt_clr'] not found — using adt.X for correlation.",
            UserWarning, stacklevel=3,
        )
        mat = adt.X
    return mat.toarray() if sp.issparse(mat) else np.asarray(mat)


def _get_rna_matrix(rna: AnnData, use_logcounts: bool) -> np.ndarray:
    """Return dense (n_cells, n_genes) log-normalised RNA matrix."""
    if use_logcounts and "logcounts" in rna.layers:
        mat = rna.layers["logcounts"]
    else:
        if use_logcounts:
            warnings.warn(
                "rna.layers['logcounts'] not found — using rna.X.",
                UserWarning, stacklevel=3,
            )
        mat = rna.X
    return mat.toarray() if sp.issparse(mat) else np.asarray(mat)


def _get_column(matrix: np.ndarray, idx: int) -> np.ndarray:
    """Return 1-D array for column idx."""
    return matrix[:, idx].ravel()
