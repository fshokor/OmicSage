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

### ✅ DEG (pipeline/modules/downstream/deg.py)
- Input: annotated AnnData with obs['cell_type_vote'] and layers['logcounts']
- Wilcoxon rank-sum via sc.tl.rank_genes_groups(), one-vs-rest per cell type
- rankby_abs=True — returns both up- and downregulated genes
- n_genes default 500 — prevents artificial cap on significant DEGs
- exclude_gene_prefixes param — post-filters RPL/RPS/MT- without biasing fold-changes
- Fallback: tries obs['cell_type_vote'], falls back to obs['leiden'] with UserWarning
- BH FDR correction; configurable min_logfc and max_pval_adj thresholds
- Provenance stored in uns['omicsage_deg']
- 11 tests passing in tests/test_deg.py

### ✅ DEG Report (reports/deg_report.py)
- max_volcano_groups default 20; truncation note + sorted by DEG count
- Direction column in Top DEGs table (▲ Up red / ▼ Down blue)
- n_genes stat card in Run Summary; exclude_prefixes info note

### ✅ GSEA (pipeline/modules/downstream/gsea.py)
- Input: deg_dict['results'] + adata (for gene universe)
- ORA via gseapy.enrichr (Fisher exact + BH correction)
- Gene sets: GO Biological Process 2023, KEGG 2021 Human, Reactome 2022
- direction param: "up" | "down" | "both"
- exclude_gene_prefixes param: filters query list only, universe unchanged
- Provenance stored in uns['omicsage_gsea']
- 8 tests passing in tests/test_gsea.py (all Enrichr calls mocked — CI-safe)

### ✅ GSEA Report (reports/gsea_report.py)
- Run summary, top pathways table, bar charts per group, bubble plot
- Direction badges: ▲ Up / ▼ Down when direction="both"

### ✅ Harmony Batch Correction (pipeline/modules/integration/harmony_correct.py) ← NEW
- Harmony integration on obs[batch_key] (default: 'batch')
- Corrected embedding stored in obsm['X_pca_harmony']
- Original UMAP preserved as obsm['X_umap_precorrection'] before overwriting
- Post-correction UMAP stored as obsm['X_umap_harmony'] (not X_umap)
- Neighbor graph recomputed on corrected embedding → uns['neighbors_harmony'],
  obsp['neighbors_harmony_connectivities'], obsp['neighbors_harmony_distances']
- Provenance stored in uns['omicsage_harmony']
  - umap_key: 'X_umap_harmony'
  - umap_precorrection_key: 'X_umap_precorrection'
- obs[batch_key] cast to str internally — safe with Categorical columns
- 12 tests passing in tests/test_harmony.py

### ✅ Harmony Report (reports/harmony_report.py) ← NEW
- Run summary: stat cards (cells, genes, batches, PCs, k, elapsed) + output key verification
- Batch composition: bar chart + table (cells per batch, % of total)
- UMAP embeddings: side-by-side raw PCA (before) vs corrected UMAP (after), coloured by batch
  + separate UMAP coloured by Harmony PC1 value (shows correction depth)
- Batch mixing metrics: per-cell same-batch neighbour fraction histogram,
  mean/median/expected stats, normalised mixing score with interpretation note
- Per-PC correction shift: bar chart of mean |X_pca − X_pca_harmony| per PC,
  top 5 most-shifted PCs table

### ✅ Notebook (notebooks/phase1_qc.ipynb)
- Step 8 (Harmony) section added — 7 cells
- Covers: load GSE194122_cite_gsea.h5ad, check batch column, run harmony_correct(),
  side-by-side PCA before vs UMAP after (using X_umap_precorrection and X_umap_harmony),
  cell-type UMAP on corrected embedding, generate HTML report, save output

## Total Tests Passing
~201 (189 pre-session + 12 new harmony tests)

## What Is NOT Built Yet
- Clustering on harmony-corrected embedding ← NEXT
- scVI batch correction (alternative to Harmony)
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
- data/processed/GSE194122_cite_gsea.h5ad
- data/processed/GSE194122_cite_harmony.h5ad   ← NEW (written by notebook Step 8)
