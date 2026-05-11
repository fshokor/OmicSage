# OmicSage — Current Status
> Last updated: 2026-05-11

## Phase
Phase 1 — Core scRNA Pipeline

## What Is Built and Tested Right Now

### ✅ Phase 0 — Foundation
- Repo structure scaffolded
- Docker base images (Python + R) — defined, not yet built locally
- CI/CD via GitHub Actions — configured
- YAML config schema — defined
- .dev_memory/ system — initialized

### ✅ Ingestion (pipeline/modules/qc/ingest.py)
- Auto-detects 10x MEX, H5, AnnData formats
- Moves normalized values out of .X, puts raw counts in .X
- Handles GSE194122 CITE-seq and multiome, GSE166635 HCC

### ✅ QC (pipeline/modules/qc/qc.py)
- Modality-aware: auto-detects GEX / ADT / ATAC from var['feature_types']
- Returns MuData with mdata["rna"], mdata["adt"], mdata["atac"] slots
- MT%, genes/cell, Scrublet doublets
- 42 tests passing in tests/test_qc.py

### ✅ Normalization (pipeline/modules/qc/normalize.py)
- Input: AnnData with raw counts in .X (from mdata["rna"])
- Saves raw counts to layers['counts']
- Normalizes to CP10K, log1p transform
- Saves log1p values to layers['logcounts'] (Seurat convention)
- HVG selection — top 2000 genes, flavor='seurat_v3' (pre-log1p on counts)
- batch_key support for per-batch HVG selection
- Provenance stored in uns['omicsage_normalization']
- 12 tests passing in tests/test_normalize.py

### ✅ Normalization Report (reports/normalization_report.py)
- Self-contained HTML report
- Figures: HVG scatter, library size violin, top 20 HVGs, gene detection rate

### ✅ Dimensionality Reduction (pipeline/modules/qc/reduce.py)
- PCA on HVG subset only (n_comps=50, svd_solver='arpack')
- Auto-select n_pcs via elbow detection (kneed) — falls back to variance threshold
- Neighbor graph (sc.pp.neighbors, n_neighbors=15), UMAP (always), t-SNE (optional)
- Provenance stored in uns['omicsage_reduce']
- 12 tests passing in tests/test_reduce.py

### ✅ Dimensionality Reduction Report (reports/reduce_report.py)
- Figures: scree plot, UMAP × QC metrics (2×2), PCA × QC, optional batch UMAP

### ✅ Leiden Clustering (pipeline/modules/clustering/cluster.py)
- Leiden sweep across configurable resolution_range
- Silhouette score per resolution; best_resolution_override for manual pinning
- obs['leiden_*'] per resolution + obs['leiden'] convenience key
- Provenance stored in uns['omicsage_cluster']
- 8 tests passing in tests/test_cluster.py

### ✅ Clustering Report (reports/cluster_report.py)
- Figures: UMAP grid per resolution, silhouette bar chart, cluster size distribution,
  optional ground-truth cell_type UMAP

### ✅ Cell-Type Annotation (pipeline/modules/annotation/annotate.py)
- Methods: CellTypist (Immune_All_High + Immune_All_Low), marker gene scoring, majority vote
- CellTypist models cached in data/references/celltypist/ (project-local)
- obs columns: celltypist_coarse, celltypist_fine, cell_type_markers,
  cell_type_groundtruth (preserved from obs['cell_type']), cell_type_vote, cell_type_confidence
- Provenance stored in uns['omicsage_annotate']
- 18 tests passing, 1 skipped in tests/test_annotate.py

### ✅ Annotation Report (reports/annotate_report.py)
- Figures: UMAP × consensus vote, UMAP × CellTypist fine, confidence distribution
- Per-cluster table with all method labels and confidence scores

### ✅ DEG (pipeline/modules/downstream/deg.py) ← UPDATED THIS SESSION
- Input: annotated AnnData with obs['cell_type_vote'] and layers['logcounts']
- Wilcoxon rank-sum via sc.tl.rank_genes_groups(), one-vs-rest per cell type
- rankby_abs=True — returns both up- and downregulated genes (critical fix)
- n_genes default raised 200→500 — prevents artificial cap on significant DEGs
- exclude_gene_prefixes param — post-filters RPL/RPS/MT- without biasing fold-changes
- Fallback: tries obs['cell_type_vote'], falls back to obs['leiden'] with UserWarning
- BH FDR correction; configurable min_logfc and max_pval_adj thresholds
- Optional pairwise mode via pairwise_groups parameter
- deg_dict: results (per-group DataFrames), summary_df (top 5 per group), provenance, pairwise
- Provenance stored in uns['omicsage_deg']
- 11 tests passing in tests/test_deg.py

### ✅ DEG Report (reports/deg_report.py) ← UPDATED THIS SESSION
- max_volcano_groups default raised 9→20 — all cell types rendered by default
- Volcano truncation note now visible in report + sorted by DEG count
- Direction column added to Top DEGs table (▲ Up red / ▼ Down blue)
- n_genes stat card added to Run Summary
- exclude_prefixes info note displayed when prefix filtering was applied

### ✅ GSEA (pipeline/modules/downstream/gsea.py) ← NEW THIS SESSION
- Input: deg_dict['results'] + adata (for gene universe)
- ORA via gseapy.enrichr (Fisher exact + BH correction)
- Gene sets: GO Biological Process 2023, KEGG 2021 Human, Reactome 2022 (configurable)
- direction param: "up" (default) | "down" | "both"
  - "up"  : upregulated query genes only
  - "down": downregulated query genes only (for suppressed pathways in cancer etc.)
  - "both": two independent ORA queries per group; results keyed as {group}__up / {group}__down
- exclude_gene_prefixes param: filters query list only, gene universe unchanged (statistically correct)
- Gene set name validation against Enrichr at runtime — warns on bad names, never crashes
- Overlap column derived from Genes column (gseapy ≥1.0 dropped Overlap natively)
- Graceful skip for groups with < min_genes DEGs — warns, never crashes
- Provenance stored in uns['omicsage_gsea']
- 8 tests passing in tests/test_gsea.py (all Enrichr calls mocked — CI-safe)

### ✅ GSEA Report (reports/gsea_report.py) ← NEW THIS SESSION
- Run summary: groups tested, direction mode, gene sets queried, total sig. pathways
- Top pathways table per group (rowspan, Genes Matched count, adj. p-value, gene list)
- Direction badges: ▲ Up (red pill) / ▼ Down (blue pill) when direction="both"
- Bar charts per group — top 10 pathways by −log₁₀(adj. p-value)
- Bubble plot — pathway × group, size = genes matched, colour = −log₁₀(adj. p-value)
  - When n_groups > max_bubble_groups: selects top N by significant pathway count (not hard skip)
  - Excluded groups listed in visible note

### ✅ Notebook (notebooks/phase1_qc.ipynb) ← UPDATED THIS SESSION
- Step 7 (GSEA) section added — 9 cells
- Covers: load GSE194122_cite_deg.h5ad, re-run deg() to recover deg_dict,
  run gsea(), sanity check (T cell activation, phagocytosis, B cell signalling),
  generate HTML report, save GSE194122_cite_gsea.h5ad

## Total Tests Passing
~189 (181 pre-session + 8 new gsea tests)

## What Is NOT Built Yet
- Harmony + scVI batch correction ← NEXT
- ScType-py + SingleR-py annotation (deferred — see docs/ANNOTATION_PLAN.md)
- Pseudobulk DEG (DESeq2-style)
- scATAC module (Phase 4)
- Spatial module (Phase 5)
- Multiome module (Phase 6)
- Streamlit UI (Phase 7)
- CLI (Phase 7)
- Quarto reports (Phase 2)

## Processed Data Files
- data/processed/GSE194122_cite_rna_qc.h5ad
- data/processed/GSE194122_cite_adt_qc.h5ad
- data/processed/GSE194122_multiome_rna_qc.h5ad
- data/processed/GSE194122_multiome_atac_qc.h5ad
- data/processed/GSE166635_HCC1_qc.h5ad
- data/processed/GSE194122_cite_normalized.h5ad
- data/processed/GSE194122_cite_reduced.h5ad
- data/processed/GSE194122_cite_clustered.h5ad
- data/processed/GSE194122_cite_annotated.h5ad
- data/processed/GSE194122_cite_deg.h5ad
- data/processed/GSE194122_cite_gsea.h5ad          ← NEW (written by notebook Step 7)
