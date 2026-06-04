"""
atac_annotate.py — ATAC gene activity scores and label transfer for OmicSage.

Two steps run in sequence:

  Step A — Gene activity scores
    Converts peak × cell count matrix into a proxy gene expression matrix by
    summing peaks that overlap each gene body + a configurable upstream promoter
    window.  Result stored in:
      obsm["gene_activity"]                 — dense float32 array (n_cells × n_genes)
      uns["gene_activity_var_names"]        — list[str] of gene names (column labels)

  Step B — Label transfer from RNA
    For multiome data, RNA and ATAC barcodes are 1:1.  The simplest and most
    robust approach is a majority-vote per Leiden cluster:
      • For each ATAC Leiden cluster, intersect its barcodes with the RNA AnnData.
      • Assign the most common obs["cell_type_vote"] label in that cluster.
      • Cells in clusters with no RNA match receive "Unknown".
    Result stored in:
      obs["atac_celltype"]                  — str label per cell

Reference
---------
  Signac GeneActivity (R):
    https://stuartlab.org/signac/articles/pbmc_vignette.html
  sc-best-practices paired integration:
    https://www.sc-best-practices.org/multimodal_integration/paired_integration.html
  Seurat v5 ATAC integration:
    https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette

Namespace conventions (do not change)
--------------------------------------
  atac.obsm["gene_activity"]             proxy gene expression (Step A)
  atac.uns["gene_activity_var_names"]    gene names for gene_activity columns
  atac.obs["atac_celltype"]             transferred label (Step B)
  atac.obs["cell_type_groundtruth"]     NeurIPS ground truth (written by atac_qc)
  atac.obs["atac_leiden"]               ATAC Leiden clusters (written by atac_reduce)

API
---
annotate_atac(
    atac,
    rna=None,
    peak_annotation=None,
    promoter_upstream_bp=2000,
    min_peaks_per_gene=1,
    leiden_key="atac_leiden",
    rna_label_key="cell_type_vote",
    inplace=False,
) -> tuple[AnnData, dict]
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData


# ---------------------------------------------------------------------------
# Step A helpers — gene activity
# ---------------------------------------------------------------------------

def _parse_peak_coords(var_names: list[str]) -> pd.DataFrame:
    """
    Parse peak var_names formatted as 'chr:start-end' or 'chr-start-end'
    into a DataFrame with columns [peak, chrom, start, end].

    Unparseable peaks are silently dropped.
    """
    records = []
    for name in var_names:
        # Support both 'chr1:100-200' and 'chr1-100-200'
        try:
            if ":" in name:
                chrom, rest = name.split(":", 1)
                start_s, end_s = rest.split("-", 1)
            else:
                parts = name.rsplit("-", 2)
                if len(parts) != 3:
                    continue
                chrom, start_s, end_s = parts
            records.append({
                "peak":  name,
                "chrom": chrom,
                "start": int(start_s),
                "end":   int(end_s),
            })
        except (ValueError, IndexError):
            continue

    if not records:
        return pd.DataFrame(columns=["peak", "chrom", "start", "end"])
    return pd.DataFrame(records)


def _build_gene_activity_from_annotation(
    atac: AnnData,
    peak_annotation: pd.DataFrame,
    promoter_upstream_bp: int,
    min_peaks_per_gene: int,
) -> tuple[np.ndarray, list[str]]:
    """
    Build gene activity matrix using a peak annotation DataFrame.

    Expected columns: peak, gene_name, distance, peak_type
    (muon format from 10x h5 files).

    A peak contributes to a gene if:
      • peak_type in {"promoter", "exon", "intron"} (gene body), OR
      • peak_type == "distal" and distance < promoter_upstream_bp, OR
      • distance == 0 (peak directly overlaps the gene)

    Falls back to all peaks in the annotation if peak_type is absent.

    Returns (matrix float32, gene_names list).
    """
    required = {"peak", "gene_name"}
    if not required.issubset(peak_annotation.columns):
        raise ValueError(
            f"peak_annotation must have columns {required}. "
            f"Found: {list(peak_annotation.columns)}"
        )

    ann = peak_annotation.copy()
    peak_var_set = set(atac.var_names)

    # Keep only peaks present in this AnnData
    ann = ann[ann["peak"].isin(peak_var_set)].copy()

    # Filter by distance/type if the columns are available
    if "peak_type" in ann.columns and "distance" in ann.columns:
        ann["distance"] = pd.to_numeric(ann["distance"], errors="coerce").fillna(0)
        body_types = {"promoter", "exon", "intron"}
        mask = (
            ann["peak_type"].isin(body_types)
            | (ann["distance"] == 0)
            | ((ann["peak_type"] == "distal") & (ann["distance"].abs() < promoter_upstream_bp))
        )
        ann = ann[mask].copy()
    elif "distance" in ann.columns:
        ann["distance"] = pd.to_numeric(ann["distance"], errors="coerce").fillna(0)
        ann = ann[ann["distance"].abs() < promoter_upstream_bp].copy()
    # else: no filter — use all peaks in annotation

    # Build peak → gene(s) mapping
    gene_to_peaks: dict[str, list[str]] = {}
    for _, row in ann.iterrows():
        gene = str(row["gene_name"])
        peak = str(row["peak"])
        gene_to_peaks.setdefault(gene, []).append(peak)

    # Apply min_peaks_per_gene filter
    gene_to_peaks = {
        g: peaks
        for g, peaks in gene_to_peaks.items()
        if len(peaks) >= min_peaks_per_gene
    }

    if not gene_to_peaks:
        return np.zeros((atac.n_obs, 0), dtype=np.float32), []

    return _sum_peaks_to_genes(atac, gene_to_peaks)


def _build_gene_activity_from_coords(
    atac: AnnData,
    promoter_upstream_bp: int,
    min_peaks_per_gene: int,
) -> tuple[np.ndarray, list[str]]:
    """
    Fallback: build gene activity using only peak coordinates.

    When peak_annotation is absent we cannot map peaks to gene names.
    Instead we create one pseudo-gene per genomic region by grouping
    adjacent peaks (within promoter_upstream_bp) on the same chromosome.
    Each group's name is 'region_{chrom}_{min_start}_{max_end}'.

    This is not biologically meaningful but produces a non-empty matrix
    that allows the pipeline to continue.  A UserWarning is emitted.
    """
    warnings.warn(
        "peak_annotation not found in atac.uns['atac']['peak_annotation'] and "
        "no peak_annotation was passed.  Falling back to coordinate-based "
        "proximity grouping.  Gene names will be synthetic region labels, "
        "not real gene names.  For production use, provide peak_annotation.",
        UserWarning,
        stacklevel=4,
    )

    peak_coords = _parse_peak_coords(list(atac.var_names))
    if peak_coords.empty:
        warnings.warn(
            "Could not parse any peak coordinates from var_names.  "
            "Gene activity matrix will be empty.",
            UserWarning,
            stacklevel=4,
        )
        return np.zeros((atac.n_obs, 0), dtype=np.float32), []

    # Simple proximity grouping: sort by chrom+start; group peaks within
    # promoter_upstream_bp of each other
    peak_coords = peak_coords.sort_values(["chrom", "start"]).reset_index(drop=True)

    gene_to_peaks: dict[str, list[str]] = {}
    current_chrom = None
    current_end   = -1
    current_group: list[str] = []
    current_start = 0

    def _flush_group() -> None:
        if not current_group:
            return
        gname = f"region_{current_chrom}_{current_start}_{current_end}"
        gene_to_peaks[gname] = list(current_group)

    for _, row in peak_coords.iterrows():
        chrom = row["chrom"]
        start = int(row["start"])
        end   = int(row["end"])
        peak  = row["peak"]

        if chrom != current_chrom or start > current_end + promoter_upstream_bp:
            _flush_group()
            current_chrom = chrom
            current_start = start
            current_end   = end
            current_group = [peak]
        else:
            current_group.append(peak)
            current_end = max(current_end, end)

    _flush_group()

    gene_to_peaks = {
        g: peaks
        for g, peaks in gene_to_peaks.items()
        if len(peaks) >= min_peaks_per_gene
    }

    if not gene_to_peaks:
        return np.zeros((atac.n_obs, 0), dtype=np.float32), []

    return _sum_peaks_to_genes(atac, gene_to_peaks)


def _sum_peaks_to_genes(
    atac: AnnData,
    gene_to_peaks: dict[str, list[str]],
) -> tuple[np.ndarray, list[str]]:
    """
    Sum peak counts across peaks assigned to each gene.

    Returns (matrix float32 shape n_cells × n_genes, sorted gene_names).
    """
    peak_index: dict[str, int] = {p: i for i, p in enumerate(atac.var_names)}

    X = atac.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    gene_names = sorted(gene_to_peaks.keys())
    out = np.zeros((atac.n_obs, len(gene_names)), dtype=np.float32)

    for col_idx, gene in enumerate(gene_names):
        peak_idxs = [peak_index[p] for p in gene_to_peaks[gene] if p in peak_index]
        if peak_idxs:
            out[:, col_idx] = X[:, peak_idxs].sum(axis=1)

    return out, gene_names


# ---------------------------------------------------------------------------
# Step B helpers — label transfer
# ---------------------------------------------------------------------------

def _majority_vote_label_transfer(
    atac: AnnData,
    rna: AnnData,
    leiden_key: str,
    rna_label_key: str,
) -> np.ndarray:
    """
    Transfer RNA cell type labels to ATAC cells via majority vote.

    Algorithm:
      1. Intersect ATAC and RNA obs_names (barcodes).
      2. For each ATAC Leiden cluster, find the majority rna_label_key
         value among cells with matching barcodes.
      3. Assign that majority label to all cells in the cluster.
      4. Cells in clusters with no barcode overlap → "Unknown".

    Returns np.ndarray of str, shape (n_atac_cells,).
    """
    if leiden_key not in atac.obs.columns:
        raise KeyError(
            f"atac.obs['{leiden_key}'] not found.  "
            "Run atac_reduce.py first to compute Leiden clusters."
        )

    if rna_label_key not in rna.obs.columns:
        raise KeyError(
            f"rna.obs['{rna_label_key}'] not found in RNA AnnData.  "
            f"Available columns: {list(rna.obs.columns)}"
        )

    # Build barcode → RNA label lookup
    rna_labels: dict[str, str] = rna.obs[rna_label_key].astype(str).to_dict()

    atac_barcodes = atac.obs_names.tolist()
    atac_leiden   = atac.obs[leiden_key].astype(str).values

    labels_out = np.full(atac.n_obs, "Unknown", dtype=object)
    unique_clusters = sorted(set(atac_leiden))

    for cluster in unique_clusters:
        cluster_mask = atac_leiden == cluster
        cluster_barcodes = [
            bc for bc, m in zip(atac_barcodes, cluster_mask) if m
        ]

        # Collect RNA labels for barcodes that exist in RNA
        rna_votes = [rna_labels[bc] for bc in cluster_barcodes if bc in rna_labels]

        if not rna_votes:
            # No barcode overlap — keep "Unknown"
            continue

        # Majority vote
        vote_counts: dict[str, int] = {}
        for label in rna_votes:
            vote_counts[label] = vote_counts.get(label, 0) + 1
        majority_label = max(vote_counts, key=lambda k: vote_counts[k])
        labels_out[cluster_mask] = majority_label

    return labels_out.astype(str)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def annotate_atac(
    atac: AnnData,
    rna: Optional[AnnData] = None,
    peak_annotation: Optional[pd.DataFrame] = None,
    promoter_upstream_bp: int = 2000,
    min_peaks_per_gene: int = 1,
    leiden_key: str = "atac_leiden",
    rna_label_key: str = "cell_type_vote",
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Compute gene activity scores and transfer RNA cell type labels to ATAC cells.

    Parameters
    ----------
    atac : AnnData
        ATAC AnnData from atac_reduce.py (must have obs[leiden_key]).
        var_names should be peak IDs in 'chr:start-end' format (10x convention).
    rna : AnnData or None
        RNA AnnData with obs[rna_label_key] (e.g. obs["cell_type_vote"]).
        When None: obs["atac_celltype"] is written as "Unknown" for all cells
        and a UserWarning is emitted.
    peak_annotation : pd.DataFrame or None
        Override for atac.uns["atac"]["peak_annotation"].
        DataFrame with columns: peak, gene_name, distance, peak_type.
        When None: falls back to uns["atac"]["peak_annotation"], then to
        coordinate-based proximity grouping.
    promoter_upstream_bp : int
        Upstream window in base pairs added to gene body for peak inclusion.
        Default: 2000.
    min_peaks_per_gene : int
        Minimum number of peaks that must overlap a gene to include it in
        the gene activity matrix.  Default: 1.
    leiden_key : str
        obs column name for ATAC Leiden cluster IDs.  Default: "atac_leiden".
    rna_label_key : str
        obs column in RNA AnnData to transfer.  Default: "cell_type_vote".
    inplace : bool
        If True, modify atac in place; default makes a copy.

    Returns
    -------
    adata : AnnData
        atac AnnData with:
          obsm["gene_activity"]           — float32 (n_cells × n_genes)
          uns["gene_activity_var_names"]  — list[str] of gene names
          obs["atac_celltype"]           — transferred label (str)
          uns["omicsage_atac_annotate"]  — provenance
    metrics : dict
        n_cells, n_peaks, n_genes_activity, atac_celltype_counts,
        n_rna_barcodes_matched, leiden_key, rna_label_key,
        promoter_upstream_bp, min_peaks_per_gene,
        peak_annotation_source, rna_provided.

    Raises
    ------
    KeyError
        If leiden_key is not found in atac.obs (Step B).
        If rna_label_key is not found in rna.obs when rna is provided.
    """
    # ------------------------------------------------------------------
    # 0. Copy or in-place
    # ------------------------------------------------------------------
    adata = atac if inplace else atac.copy()

    # ------------------------------------------------------------------
    # A. Resolve peak annotation source
    # ------------------------------------------------------------------
    ann_df: Optional[pd.DataFrame] = None
    peak_annotation_source: str = "none"

    if peak_annotation is not None:
        ann_df = peak_annotation
        peak_annotation_source = "provided"
    else:
        # Try to pull from uns["atac"]["peak_annotation"]
        atac_uns = adata.uns.get("atac", {})
        if isinstance(atac_uns, dict) and "peak_annotation" in atac_uns:
            candidate = atac_uns["peak_annotation"]
            if isinstance(candidate, pd.DataFrame) and not candidate.empty:
                ann_df = candidate
                peak_annotation_source = "uns"

    # ------------------------------------------------------------------
    # A1. Compute gene activity matrix
    # ------------------------------------------------------------------
    if ann_df is not None:
        gene_activity_matrix, gene_names = _build_gene_activity_from_annotation(
            adata,
            ann_df,
            promoter_upstream_bp=promoter_upstream_bp,
            min_peaks_per_gene=min_peaks_per_gene,
        )
    else:
        gene_activity_matrix, gene_names = _build_gene_activity_from_coords(
            adata,
            promoter_upstream_bp=promoter_upstream_bp,
            min_peaks_per_gene=min_peaks_per_gene,
        )
        if peak_annotation_source == "none":
            peak_annotation_source = "coordinate_fallback"

    adata.obsm["gene_activity"]          = gene_activity_matrix
    adata.uns["gene_activity_var_names"] = gene_names
    n_genes_activity                     = len(gene_names)

    # ------------------------------------------------------------------
    # B. Label transfer from RNA
    # ------------------------------------------------------------------
    n_rna_barcodes_matched = 0

    if rna is None:
        warnings.warn(
            "No RNA AnnData provided.  "
            "obs['atac_celltype'] will be 'Unknown' for all cells.  "
            "Pass rna= to enable label transfer.",
            UserWarning,
            stacklevel=2,
        )
        adata.obs["atac_celltype"] = "Unknown"
    else:
        transferred = _majority_vote_label_transfer(
            adata,
            rna,
            leiden_key=leiden_key,
            rna_label_key=rna_label_key,
        )
        adata.obs["atac_celltype"] = transferred

        # Count matched barcodes for metrics
        rna_barcodes = set(rna.obs_names)
        n_rna_barcodes_matched = int(
            sum(1 for bc in adata.obs_names if bc in rna_barcodes)
        )

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------
    adata.uns["omicsage_atac_annotate"] = {
        "module":    "atac_annotate",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "promoter_upstream_bp":   promoter_upstream_bp,
            "min_peaks_per_gene":     min_peaks_per_gene,
            "leiden_key":             leiden_key,
            "rna_label_key":          rna_label_key,
            "peak_annotation_source": peak_annotation_source,
            "rna_provided":           rna is not None,
        },
        "outputs": {
            "gene_activity_key":      "gene_activity",
            "gene_activity_var_names_key": "gene_activity_var_names",
            "n_genes_activity":       n_genes_activity,
            "celltype_key":           "atac_celltype",
            "n_rna_barcodes_matched": n_rna_barcodes_matched,
        },
    }

    # ------------------------------------------------------------------
    # Metrics
    # ------------------------------------------------------------------
    celltype_counts: dict[str, int] = (
        adata.obs["atac_celltype"]
        .value_counts()
        .to_dict()
    )
    # Ensure JSON-serialisable int values
    celltype_counts = {str(k): int(v) for k, v in celltype_counts.items()}

    metrics: dict = {
        "n_cells":                int(adata.n_obs),
        "n_peaks":                int(adata.n_vars),
        "n_genes_activity":       n_genes_activity,
        "atac_celltype_counts":   celltype_counts,
        "n_rna_barcodes_matched": n_rna_barcodes_matched,
        "leiden_key":             leiden_key,
        "rna_label_key":          rna_label_key,
        "promoter_upstream_bp":   promoter_upstream_bp,
        "min_peaks_per_gene":     min_peaks_per_gene,
        "peak_annotation_source": peak_annotation_source,
        "rna_provided":           rna is not None,
    }

    return adata, metrics
