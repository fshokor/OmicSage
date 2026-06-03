"""
cite_corr.py — Protein-RNA Correlation module for OmicSage
Phase 4, Step 9

Computes per-protein Spearman correlation between ADT expression and the
corresponding gene's RNA expression across all cells.

For each protein in the ADT panel, the module looks for a gene in the RNA
layer with the same name (case-insensitive, after stripping common suffixes
like "-TotalSeqA"). Matched pairs are tested with scipy.stats.spearmanr.

Biological interpretation
-------------------------
High correlation (r > 0.4):
    Protein abundance is tightly coupled to mRNA — standard housekeeping
    or lineage markers (e.g. CD3E, CD19, CD14).

Low or zero correlation (|r| < 0.2):
    Post-transcriptional regulation — protein abundance is controlled at
    the translation, secretion, or membrane-trafficking level rather than
    transcription. Biologically interesting; may indicate regulatory events.

Negative correlation:
    Rare but meaningful. Can indicate feedback inhibition, splice variants,
    or technical artefact (cross-reactive antibody binding a different epitope).

Low-correlation proteins are often the most scientifically interesting because
they highlight biology that neither modality captures alone.

Input
-----
MuData from cite_06 integration step (preferred).
Requires:
  mdata["adt"].layers["adt_clr"]   — CLR-normalised ADT (float matrix)
  mdata["rna"].X or layers["logcounts"]  — log-normalised RNA
  mdata["adt"].var_names            — protein names (e.g. "CD3E", "CD19")
  mdata["rna"].var_names            — gene symbols (e.g. "CD3E", "CD19")

AnnData fallback: if a bare AnnData (cite_05) is passed, the module raises
a ValueError immediately — correlation requires both modalities.

Matching strategy
-----------------
1. Exact match (case-insensitive): "CD3E" → "CD3E"
2. Suffix stripping: "CD3E-TotalSeqA" → "CD3E" → match if "CD3E" in RNA
3. No match: protein skipped, logged in corr_dict["unmatched"]

Only proteins with >= min_cells cells expressing both modalities
(ADT CLR > 0 AND RNA > 0) are tested.

Outputs written to data.uns
---------------------------
uns["omicsage_cite_corr"]   Provenance block

Returns
-------
(MuData, corr_dict)

corr_dict keys
--------------
  "results"    : DataFrame(protein, gene, r, pval, pval_adj, n_cells, matched_by)
  "matched"    : list[tuple[str, str]] — (protein, gene) pairs tested
  "unmatched"  : list[str] — protein names with no RNA match
  "provenance" : same as uns["omicsage_cite_corr"]

Usage
-----
    from pipeline.modules.cite.cite_corr import cite_corr

    mdata_corr, corr_dict = cite_corr(
        mdata,
        method="spearman",
        min_cells=30,
        matched_only=True,
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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cite_corr(
    data: Union["MuData", AnnData],
    method: str = "spearman",
    min_cells: int = 30,
    use_logcounts: bool = True,
    inplace: bool = False,
) -> tuple[Union["MuData", AnnData], dict]:
    """
    Compute per-protein Spearman correlation between ADT and RNA expression.

    Parameters
    ----------
    data : MuData
        MuData from cite_06 integration step.
        Requires mdata["adt"] and mdata["rna"] modalities.
        AnnData input raises ValueError immediately.
    method : str
        Correlation method. Currently only "spearman" is supported.
        "pearson" support is planned but not implemented — CLR values are
        not normally distributed so Pearson is inappropriate without further
        transformation.
    min_cells : int
        Minimum number of cells expressing both the protein (ADT CLR > 0)
        and the gene (RNA > 0) for a pair to be tested.
        Pairs below this threshold are skipped and noted in corr_dict["skipped"].
        Default: 30.
    use_logcounts : bool
        If True (default), uses mdata["rna"].layers["logcounts"] for RNA
        expression. Falls back to mdata["rna"].X with a UserWarning if the
        layer is absent.
    inplace : bool
        If False (default), operates on a copy.

    Returns
    -------
    data : MuData
        Input MuData with uns["omicsage_cite_corr"] populated.
    corr_dict : dict
        Keys: results, matched, unmatched, skipped, provenance.
        See module docstring for full description.

    Raises
    ------
    ValueError
        If data is an AnnData (both modalities required).
    ValueError
        If method is not "spearman".
    """
    # ------------------------------------------------------------------
    # 0. Validate input type
    # ------------------------------------------------------------------
    if not (_MUDATA_AVAILABLE and isinstance(data, MuData)):
        raise ValueError(
            "cite_corr requires a MuData object with both 'rna' and 'adt' "
            "modalities (from cite_06 integration step). "
            "AnnData (single-modality) input is not supported — "
            "protein-RNA correlation requires both layers."
        )

    if method != "spearman":
        raise ValueError(
            f"method='{method}' is not supported. "
            "Only 'spearman' is currently implemented. "
            "CLR-normalised ADT values are not normally distributed, "
            "making Pearson correlation inappropriate without further "
            "transformation."
        )

    if not inplace:
        data = data.copy()

    rna = data["rna"]
    adt = data["adt"]

    # ------------------------------------------------------------------
    # 1. Extract expression matrices
    # ------------------------------------------------------------------
    adt_matrix = _get_adt_matrix(adt)        # (n_cells, n_proteins)
    rna_matrix = _get_rna_matrix(rna, use_logcounts)  # (n_cells, n_genes)

    adt_proteins = list(adt.var_names)
    rna_genes    = list(rna.var_names)

    # Build gene lookup: lowercase → original index for fast matching
    rna_gene_lower = {g.lower(): i for i, g in enumerate(rna_genes)}

    # ------------------------------------------------------------------
    # 2. Match proteins to genes
    # ------------------------------------------------------------------
    matched:   list[tuple[str, str, int, int]] = []  # (protein, gene, prot_idx, gene_idx)
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
            "Check that mdata['adt'].var_names use standard gene-symbol names "
            "(e.g. 'CD3E', 'CD19') or names that can be stripped to gene symbols.",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 3. Compute correlations
    # ------------------------------------------------------------------
    rows: list[dict] = []
    skipped: list[dict] = []

    for protein, gene, prot_idx, gene_idx, matched_by in matched:
        x_prot = _get_column(adt_matrix, prot_idx)  # ADT CLR vector
        x_rna  = _get_column(rna_matrix, gene_idx)  # RNA log vector

        # Filter to cells expressing both
        mask = (x_prot > 0) & (x_rna > 0)
        n_cells_expr = int(mask.sum())

        if n_cells_expr < min_cells:
            skipped.append({
                "protein":  protein,
                "gene":     gene,
                "n_cells":  n_cells_expr,
                "reason":   f"< min_cells ({min_cells})",
            })
            continue

        r, pval = stats.spearmanr(x_prot[mask], x_rna[mask])

        rows.append({
            "protein":    protein,
            "gene":       gene,
            "r":          float(r),
            "pval":       float(pval),
            "n_cells":    n_cells_expr,
            "matched_by": matched_by,
        })

    # ------------------------------------------------------------------
    # 4. Multiple testing correction (BH FDR)
    # ------------------------------------------------------------------
    if rows:
        results_df = pd.DataFrame(rows)
        if _STATSMODELS_AVAILABLE and len(results_df) > 1:
            _, pval_adj, _, _ = multipletests(
                results_df["pval"].values,
                method="fdr_bh",
            )
            results_df["pval_adj"] = pval_adj
        else:
            # Fallback: no correction, copy raw p-values
            results_df["pval_adj"] = results_df["pval"]
            if not _STATSMODELS_AVAILABLE:
                warnings.warn(
                    "statsmodels not available — BH FDR correction skipped. "
                    "Install with: pip install statsmodels",
                    UserWarning,
                    stacklevel=2,
                )
        results_df = results_df.sort_values("r", ascending=False).reset_index(drop=True)
    else:
        results_df = pd.DataFrame(
            columns=["protein", "gene", "r", "pval", "pval_adj",
                     "n_cells", "matched_by"]
        )

    # ------------------------------------------------------------------
    # 5. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "module":           "cite_corr",
        "timestamp":        datetime.now(timezone.utc).isoformat(),
        "method":           method,
        "min_cells":        min_cells,
        "use_logcounts":    use_logcounts,
        "n_proteins_panel": len(adt_proteins),
        "n_matched":        len(matched),
        "n_unmatched":      len(unmatched),
        "n_tested":         len(rows),
        "n_skipped":        len(skipped),
        "bh_correction":    _STATSMODELS_AVAILABLE and len(rows) > 1,
        "unmatched_proteins": unmatched,
    }
    data.uns["omicsage_cite_corr"] = provenance

    # Persist results DataFrame in uns as a JSON string so it survives
    # the h5mu round-trip. A plain records list fails because h5py cannot
    # serialise mixed-type Python objects (float, int, str, NaN) when
    # anndata writes the list as a vlen string array. A single JSON string
    # is always safe — anndata writes it as a scalar string dataset.
    if not results_df.empty:
        data.uns["omicsage_cite_corr_results"] = results_df.to_json(orient="records")
    else:
        data.uns["omicsage_cite_corr_results"] = "[]"

    # ------------------------------------------------------------------
    # 6. Assemble return dict
    # ------------------------------------------------------------------
    corr_dict = {
        "results":    results_df,
        "matched":    [(p, g) for p, g, _, _, _ in matched],
        "unmatched":  unmatched,
        "skipped":    skipped,
        "provenance": provenance,
    }

    return data, corr_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

# Common suffixes that vendors append to protein names in ADT panels.
# Stripping these often recovers the underlying gene symbol.
_STRIP_SUFFIXES = re.compile(
    r"[-_](TotalSeqA|TotalSeqB|TotalSeqC|TotalSeq|"
    r"BV421|BV711|APC|PE|FITC|BioLegend|[0-9]+).*$",
    re.IGNORECASE,
)


def _match_protein_to_gene(
    protein: str,
    rna_gene_lower: dict[str, int],
    rna_genes: list[str],
) -> tuple[Optional[int], Optional[str]]:
    """
    Try to match a protein name to a gene symbol.

    Returns (gene_index, matched_by) or (None, None).

    matched_by values:
      "exact"   — protein name matched gene symbol directly (case-insensitive)
      "stripped" — match found after removing vendor suffixes
    """
    # Step 1: exact match (case-insensitive)
    key = protein.lower()
    if key in rna_gene_lower:
        return rna_gene_lower[key], "exact"

    # Step 2: strip vendor suffixes and retry
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
            "adt.layers['adt_clr'] not found — using adt.X for correlation. "
            "Run adt_normalize first for best results.",
            UserWarning,
            stacklevel=3,
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
                "rna.layers['logcounts'] not found — using rna.X for correlation. "
                "Ensure rna.X contains log-normalised values.",
                UserWarning,
                stacklevel=3,
            )
        mat = rna.X
    return mat.toarray() if sp.issparse(mat) else np.asarray(mat)


def _get_column(matrix: np.ndarray, idx: int) -> np.ndarray:
    """Return 1-D array for column idx, handling both dense and row-major."""
    return matrix[:, idx].ravel()
