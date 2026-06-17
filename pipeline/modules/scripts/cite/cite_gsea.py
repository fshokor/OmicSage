"""
cite_gsea.py — CITE-seq GSEA module for OmicSage
Phase 4, Step 8

Runs over-representation analysis (ORA) on the cross-modal RNA DEG results
produced by cite_deg (Step 7). Uses the same gseapy-backed gsea() function
as the RNA pipeline — no new statistical logic, just a CITE-seq-aware wrapper.

Why cross-modal only
--------------------
GSEA requires gene symbols that map to GO/KEGG/Reactome gene sets. The ADT
layer contains protein names (CD3, CD19, ...) not gene symbols — enrichment
databases don't cover these. The correct approach for CITE-seq is:

  ADT clusters define immunophenotypic cell populations
      ↓
  RNA DEG on those populations (cite_deg cross-modal output)
      ↓
  GSEA on the RNA gene lists (this module)

This gives pathway enrichment for transcriptionally-defined programs of
surface-phenotype-pure populations — cleaner than RNA-cluster-based GSEA
because ADT defines cell types with less ambiguity than transcriptomics alone.

Input
-----
  data         : MuData from cite_deg (uns["omicsage_cite_deg"] must exist)
  cite_deg_dict: dict returned by cite_deg() — must contain "rna_crossmodal"

Output
------
  data.uns["omicsage_cite_gsea"]  — provenance block
  cite_gsea_dict keys:
    "results"     — {cell_type: DataFrame} — pathway results per cell type
    "summary_df"  — top 3 pathways per cell type
    "provenance"  — same as uns block
    "skipped"     — groups skipped (< min_genes DEGs)

Usage
-----
    from pipeline.modules.scripts.cite.cite_gsea import cite_gsea

    mdata_gsea, gsea_dict = cite_gsea(
        mdata,
        cite_deg_dict=deg_dict,
    )
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional, Union

import pandas as pd
from anndata import AnnData

try:
    from mudata import MuData
    _MUDATA_AVAILABLE = True
except ImportError:
    _MUDATA_AVAILABLE = False
    MuData = None  # type: ignore

# Re-use the existing RNA pipeline GSEA function
from pipeline.modules.scripts.downstream.gsea import gsea as _rna_gsea


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cite_gsea(
    data: Union["MuData", AnnData],
    cite_deg_dict: dict,
    gene_sets: Optional[list[str]] = None,
    min_logfc: float = 0.25,
    max_pval_adj: float = 0.05,
    top_n_genes: Optional[int] = None,
    min_genes: int = 5,
    organism: str = "human",
    direction: str = "up",
    exclude_gene_prefixes: Optional[list[str]] = None,
    inplace: bool = False,
) -> tuple[Union["MuData", AnnData], dict]:
    """
    Run ORA pathway enrichment on cross-modal RNA DEG results from cite_deg.

    Wraps the existing RNA pipeline gsea() function, adapting its inputs and
    outputs for the CITE-seq MuData context.

    Parameters
    ----------
    data : MuData or AnnData
        Object returned by cite_deg (Step 7). Must contain
        uns["omicsage_cite_deg"].
    cite_deg_dict : dict
        Dict returned by cite_deg(). Must contain "rna_crossmodal" key with
        {cell_type: DataFrame(gene, score, pval, logfc, pval_adj)}.
        If "rna_crossmodal" is empty (AnnData fallback path), raises
        ValueError with a clear message.
    gene_sets : list of str, optional
        Enrichr gene set library names. Defaults to:
        ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"].
    min_logfc : float
        Minimum absolute log2FC for a gene to enter the ORA query list.
    max_pval_adj : float
        Maximum BH-adjusted p-value for a gene to enter the ORA query list.
    top_n_genes : int or None
        Use only the top N genes per group. None = all passing filters.
    min_genes : int
        Minimum query genes required to run enrichment for a group.
        Groups below this threshold are skipped with a UserWarning.
    organism : str
        Enrichr organism string. Default: "human".
    direction : str
        "up" | "down" | "both". Default: "up".
        See gsea() docstring for full description.
    exclude_gene_prefixes : list of str, optional
        Gene prefixes to exclude from ORA query list.
        Example: ["RPL", "RPS", "MT-"].
    inplace : bool
        If False (default), operates on a copy.

    Returns
    -------
    data : MuData or AnnData
        Input object with uns["omicsage_cite_gsea"] populated.
    cite_gsea_dict : dict
        See module docstring for key descriptions.

    Raises
    ------
    ValueError
        If cite_deg_dict["rna_crossmodal"] is empty (cite_deg ran in
        AnnData fallback mode — cross-modal DEG was not computed).
    """
    is_mudata = _MUDATA_AVAILABLE and isinstance(data, MuData)

    if not inplace:
        data = data.copy()

    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    rna_results = cite_deg_dict.get("rna_crossmodal", {})
    if not rna_results:
        raise ValueError(
            "cite_gsea requires cross-modal RNA DEG results from cite_deg. "
            "cite_deg_dict['rna_crossmodal'] is empty — this happens when "
            "cite_deg was run on an AnnData (cite_05) instead of a MuData "
            "(cite_06). Re-run cite_deg with the cite_06 MuData checkpoint "
            "to enable cross-modal DEG and GSEA."
        )

    if "omicsage_cite_deg" not in (
        data.uns if not is_mudata else data.uns
    ):
        warnings.warn(
            "uns['omicsage_cite_deg'] not found. "
            "Ensure cite_deg was run before cite_gsea.",
            UserWarning,
            stacklevel=2,
        )

    # ------------------------------------------------------------------
    # 2. Build a minimal AnnData stub for the RNA gene universe
    #    gsea() needs adata.var_names for background gene set
    # ------------------------------------------------------------------
    rna_adata = _get_rna_adata(data, is_mudata)

    # ------------------------------------------------------------------
    # 3. Build deg_dict in the format gsea() expects
    #    gsea() needs: {"results": {group: DataFrame}, "provenance": {...}}
    # ------------------------------------------------------------------
    deg_provenance = cite_deg_dict.get("provenance", {})
    adapted_deg_dict = {
        "results": rna_results,
        "provenance": {
            "groupby": deg_provenance.get("groupby", "adt_celltype_manual"),
            **deg_provenance,
        },
    }

    # ------------------------------------------------------------------
    # 4. Run ORA via existing RNA gsea() function
    # ------------------------------------------------------------------
    _, gsea_result = _rna_gsea(
        rna_adata,
        deg_dict=adapted_deg_dict,
        gene_sets=gene_sets,
        min_logfc=min_logfc,
        max_pval_adj=max_pval_adj,
        top_n_genes=top_n_genes,
        min_genes=min_genes,
        organism=organism,
        direction=direction,
        exclude_gene_prefixes=exclude_gene_prefixes,
        inplace=False,
    )

    # ------------------------------------------------------------------
    # 5. Provenance — tag as cite_gsea, not gsea
    # ------------------------------------------------------------------
    provenance = {
        **gsea_result["provenance"],
        "module": "cite_gsea",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_type": cite_deg_dict.get("input_type", "unknown"),
        "source_deg_module": "cite_deg",
        "note": (
            "ORA run on cross-modal RNA DEG gene lists "
            "(ADT-defined cell types → RNA transcriptional programs)"
        ),
    }

    # Write to top-level uns
    data.uns["omicsage_cite_gsea"] = provenance

    # ------------------------------------------------------------------
    # 6. Assemble return dict
    # ------------------------------------------------------------------
    cite_gsea_dict = {
        "results": gsea_result["results"],
        "summary_df": gsea_result["summary_df"],
        "provenance": provenance,
        "skipped": gsea_result.get("skipped", []),
    }

    return data, cite_gsea_dict


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _get_rna_adata(
    data: Union["MuData", AnnData],
    is_mudata: bool,
) -> AnnData:
    """
    Return an AnnData with var_names representing the RNA gene universe.

    gsea() uses adata.var_names as the ORA background gene set. We pass
    the RNA modality from the MuData, or the object itself if AnnData.
    """
    if is_mudata:
        return data["rna"]

    # AnnData fallback — shouldn't reach here (we raise earlier) but be safe
    warnings.warn(
        "cite_gsea: using AnnData var_names as gene universe. "
        "Ensure this is the RNA AnnData, not the ADT AnnData.",
        UserWarning,
        stacklevel=3,
    )
    return data
