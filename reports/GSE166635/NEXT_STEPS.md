# OmicSage — Suggested Next Steps

*Generated: 2026-05-16 17:56 UTC*


## Overview

The prioritization logic is based on the biological question and the coherence review flags. The most impactful suggestion is to re-cluster and re-annotate cluster DC2, as this will provide accurate information about the cell type identity of cells in this cluster. The subsequent suggestions focus on characterizing the immune and stromal microenvironment, identifying immunosuppressive cell populations, and determining how the tumour microenvironment differs from adjacent normal tissue.

## 1. Cell-cell communication analysis

**Rationale**: Immune and non-immune cell populations co-exist in the dataset. Ligand-receptor analysis will reveal signalling interactions between the microenvironment compartments.

**Expected output**: Ranked ligand-receptor pairs between cell type pairs, with pathway-level summaries of dominant signalling axes.

**Relevant tool**: CellChat / LIANA


## 2. Sub-clustering of heterogeneous populations

**Rationale**: Coherence review identified cluster(s) 10, 11, 12 as containing mixed cell type markers, suggesting unresolved heterogeneity that finer resolution clustering can reveal.

**Expected output**: Sub-clusters within the flagged populations, with marker genes distinguishing each sub-population.

**Relevant tool**: Leiden clustering at higher resolution (scanpy.tl.leiden)


## 3. Pseudobulk differential expression

**Rationale**: The dataset has 2 conditions. Pseudobulk DEG aggregates counts per donor before testing, providing better control of donor-level variance than single-cell Wilcoxon tests.

**Expected output**: Per-condition DEG tables with effect sizes and adjusted p-values, suitable for downstream pathway enrichment.

**Relevant tool**: PyDESeq2 / DESeq2 (R)


## 4. Re-cluster and re-annotate cluster DC2

**Rationale**: The coherence review flags a warning about the inconsistent upregulation of LST1, CST3, and COTL1 genes in cluster DC2. Verifying the cell type assignment for this cluster is crucial to understand its role in the tumour microenvironment.

**Expected output**: A revised cluster assignment that accurately reflects the biological identity of cells in cluster DC2

**Relevant tool**: Seurat or Monocle


## 5. Compare immune and stromal cell populations between tumour and adjacent normal tissue

**Rationale**: The biological question focuses on characterizing the immune and stromal microenvironment of hepatocellular carcinoma. Comparing these populations between tumour and adjacent normal tissue will provide valuable insights into how the tumour microenvironment differs at single-cell resolution.

**Expected output**: A comprehensive comparison of immune and stromal cell populations between tumour and adjacent normal tissue

**Relevant tool**: Seurat or Monocle


## 6. Identify immunosuppressive cell populations using marker genes

**Rationale**: The biological question aims to identify immunosuppressive cell populations. Using marker genes, such as PD-1 or CTLA4, can help pinpoint these cells and their potential therapeutic targets.

**Expected output**: A list of immunosuppressive cell populations and their corresponding marker genes

**Relevant tool**: Seurat or Monocle


## 7. Analyze the expression of immune checkpoint molecules on immune cells

**Rationale**: The biological question focuses on characterizing the immune microenvironment. Analyzing the expression of immune checkpoint molecules, such as PD-1 and CTLA4, on immune cells can provide insights into their potential role in immunosuppression.

**Expected output**: The expression profile of immune checkpoint molecules on immune cells

**Relevant tool**: Seurat or Monocle


## 8. Compare the expression of genes involved in immune regulation between tumour and adjacent normal tissue

**Rationale**: The biological question aims to determine how the tumour microenvironment differs from adjacent normal tissue. Comparing the expression of genes involved in immune regulation can provide valuable insights into these differences.

**Expected output**: A comparison of gene expression involved in immune regulation between tumour and adjacent normal tissue

**Relevant tool**: Seurat or Monocle


## 9. Identify potential therapeutic targets using single-cell data

**Rationale**: The biological question aims to identify potential therapeutic targets. Using single-cell data, researchers can identify specific cell populations and their corresponding marker genes that may be targeted therapeutically.

**Expected output**: A list of potential therapeutic targets identified from single-cell data

**Relevant tool**: Seurat or Monocle

