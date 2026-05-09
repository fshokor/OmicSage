# OmicSage — Current Status
> Last updated: 2026-05-09

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

### ✅ DEG (pipeline/modules/qc/deg.py)  ← NEW THIS SESSION
- Input: annotated AnnData with obs['cell_type_vote'] and layers['logcounts']
- Wilcoxon rank-sum via sc.tl.rank_genes_groups(), one-vs-rest per cell type
- Fallback: tries obs['cell_type_vote'], falls back to obs['leiden'] with UserWarning
- BH FDR correction; configurable min_logfc and max_pval_adj thresholds
- Threshold filtering applied post-extraction (uns['rank_genes_groups'] always has full n_genes)
- pts=True passed to scanpy — fraction expressing stored at no cost
- Optional pairwise mode via pairwise_groups parameter
- Warns for groups with < 10 cells (unreliable Wilcoxon)
- deg_dict: results (per-group DataFrames), summary_df (top 5 per group), provenance, pairwise
- Provenance stored in uns['omicsage_deg']
- inplace=False default
- 11 tests passing in tests/test_deg.py

### ✅ DEG Report (reports/deg_report.py)  ← NEW THIS SESSION
- Self-contained HTML report
- Sections: run summary cards + per-group DEG counts, top DEGs table (rowspan, logFC coloured),
  volcano plots (one per group, capped at max_volcano_groups=9), dot plot (sc.pl.dotplot)
- All sections fail gracefully — one broken plot never kills the report

### ✅ Notebook (notebooks/phase1_qc.ipynb)  ← UPDATED THIS SESSION
- Step 6 (DEG) section added — 9 cells
- Covers: load annotated AnnData, run deg(), sanity checks, biological marker check
  (CD3D → T cell, CD14 → Monocyte, MS4A1 → B cell), report generation, save
- Output: data/processed/GSE194122_cite_deg.h5ad

### ✅ MODULE_DOCS.md  ← UPDATED THIS SESSION
- Modules 11–14 documented: annotate.py, annotate_report.py, deg.py, deg_report.py
- Data flow diagram extended to DEG → GSEA
- Tests table updated with accurate counts

## Total Tests Passing
~181 (170 pre-session + 11 new deg tests)

## What Is NOT Built Yet
- GSEA module (blocked until deg.py tests confirmed passing on real data)
- ScType-py + SingleR-py annotation (deferred — see docs/ANNOTATION_PLAN.md)
- Harmony + scVI batch correction
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
- data/processed/GSE194122_cite_deg.h5ad          ← NEW (written by notebook Step 6)
