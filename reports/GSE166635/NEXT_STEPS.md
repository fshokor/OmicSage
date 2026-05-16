# OmicSage — Suggested Next Steps

*Generated: 2026-05-16 12:37 UTC*


## Overview

Rule-based analysis identified potential downstream steps based on detected cell types, experimental design, and coherence review flags.

## 1. Trajectory analysis

**Rationale**: Both progenitor and mature cell types are present in the dataset, suggesting a differentiation continuum that trajectory analysis can resolve.

**Expected output**: Pseudotime ordering of cells along differentiation axes, with branch points identifying fate decisions.

**Relevant tool**: Slingshot / PAGA (scVelo for RNA velocity if spliced counts available)


## 2. Cell-cell communication analysis

**Rationale**: Immune and non-immune cell populations co-exist in the dataset. Ligand-receptor analysis will reveal signalling interactions between the microenvironment compartments.

**Expected output**: Ranked ligand-receptor pairs between cell type pairs, with pathway-level summaries of dominant signalling axes.

**Relevant tool**: CellChat / LIANA


## 3. Sub-clustering of heterogeneous populations

**Rationale**: Coherence review identified cluster(s) 0, 10 as containing mixed cell type markers, suggesting unresolved heterogeneity that finer resolution clustering can reveal.

**Expected output**: Sub-clusters within the flagged populations, with marker genes distinguishing each sub-population.

**Relevant tool**: Leiden clustering at higher resolution (scanpy.tl.leiden)


## 4. Pseudobulk differential expression

**Rationale**: The dataset has 2 conditions. Pseudobulk DEG aggregates counts per donor before testing, providing better control of donor-level variance than single-cell Wilcoxon tests.

**Expected output**: Per-condition DEG tables with effect sizes and adjusted p-values, suitable for downstream pathway enrichment.

**Relevant tool**: PyDESeq2 / DESeq2 (R)

