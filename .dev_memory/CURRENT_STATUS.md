# OmicSage — Current Status
> Last updated: 2026-05-12 (session 7)

## Phase
Phase 2 — Report Engine (Phase 1 complete ✅)

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
- NEW: load_dataset_dir() — scans parent folder, loads all MTX subfolders,
  concatenates with obs['sample'] + obs['batch'] per subfolder name
- NEW: load_dataset() auto-routes to load_dataset_dir() when given a parent
  folder containing MTX subfolders (no config change needed)

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
- neighbors_key param — routes to Harmony graph when set
- cluster_key param — stores results in obs['leiden_harmony']
- compute_ari(adata, key_a, key_b) — ARI comparison between any two obs columns
- Provenance stored in uns['omicsage_cluster']
- 16 tests passing in tests/test_cluster.py

### ✅ Clustering Report (reports/cluster_report.py)
- Figures: UMAP grid per resolution, silhouette bar chart, cluster size distribution,
  optional ground-truth cell_type UMAP

### ✅ Cell-Type Annotation (pipeline/modules/annotation/annotate.py)
- Methods: CellTypist (Immune_All_High + Immune_All_Low), marker gene scoring, majority vote
- CellTypist models cached in data/references/celltypist/
- obs columns: celltypist_coarse, celltypist_fine, cell_type_markers,
  cell_type_groundtruth, cell_type_vote, cell_type_confidence
- Provenance stored in uns['omicsage_annotate']
- 18 tests passing, 1 skipped in tests/test_annotate.py

### ✅ Annotation Report (reports/annotate_report.py)
- Figures: UMAP × consensus vote, UMAP × CellTypist fine, confidence distribution
- Per-cluster table with all method labels and confidence scores

### ✅ DEG (pipeline/modules/downstream/deg.py)
- Wilcoxon rank-sum, one-vs-rest, rankby_abs=True, BH correction
- n_genes default 500, exclude_gene_prefixes param
- Provenance stored in uns['omicsage_deg']
- 11 tests passing in tests/test_deg.py

### ✅ DEG Report (reports/deg_report.py)
- Volcano plots + dot plot + summary table + direction column

### ✅ GSEA (pipeline/modules/downstream/gsea.py)
- ORA via gseapy.enrichr — GO BP / KEGG / Reactome
- direction param: "up" | "down" | "both"
- Provenance stored in uns['omicsage_gsea']
- 8 tests passing in tests/test_gsea.py (Enrichr calls mocked — CI-safe)

### ✅ GSEA Report (reports/gsea_report.py)
- Bar charts + bubble plot + direction badges

### ✅ Harmony Batch Correction (pipeline/modules/integration/harmony_correct.py)
- obsm['X_pca_harmony'], obsm['X_umap_precorrection'], obsm['X_umap_harmony']
- neighbors_harmony graph in uns + obsp
- Provenance stored in uns['omicsage_harmony']
- 13 tests passing in tests/test_harmony.py

### ✅ Harmony Report (reports/harmony_report.py)
- Batch composition, mixing metrics, PC shift, UMAP comparison

### ✅ Pseudobulk DEG (pipeline/modules/downstream/pseudobulk_deg.py)
- DESeq2 Wald tests via pydeseq2, one-vs-rest per cell type
- Aggregate layers['counts'] per (cell_type, donor)
- min_cells + min_samples filters with graceful skip + UserWarning
- Output schema identical to deg.py deg_dict
- Provenance stored in uns['omicsage_pseudobulk_deg']
- 14 tests passing in tests/test_pseudobulk_deg.py

### ✅ Pseudobulk DEG Report (reports/pseudobulk_deg_report.py)
- Skipped groups section + pseudobulk-specific stat cards

### ✅ Notebook (notebooks/phase1_qc.ipynb)
- Steps 1–10 complete

### ✅ MILESTONE: Wang et al. 2025 HCC Benchmark ← DONE THIS SESSION
- Full Phase 1 pipeline run on GSE166635 (HCC1 normal + HCC2 tumour)
- Cell types identified: Hepatocytes, T cells, Macrophages, Endothelial,
  Fibroblasts, B cells
- Known HCC markers recovered in DEG results (AFP, GPC3, EPCAM in hepatocytes;
  CD3D in T cells; CD68 in macrophages)
- Liver metabolism and immune pathway terms in GSEA results
- Reports generated: reports/GSE166635/ (all 9 steps, pseudobulk skipped)
- Processed files: data/processed/GSE166635/ (steps 01–09)

### ✅ Generic Pipeline Runner (run_pipeline.py)
- Config-driven, --from-step / --to-step / --step flags
- Validation at startup, caching, resolution_override support
- Fresh-run validation bug fixed

### ✅ Config System (config/)
- config/schema.yaml — platform master schema
- config/runs/GSE194122.yaml, GSE166635.yaml, GSE194122_multiome.yaml

## Phase 2 — What Is Being Built Next
- Quarto report templates (HTML + PDF) replacing hand-coded HTML reports
- python-pptx slide deck generator
- Auto-methods text from provenance metadata
- Goal: one command → complete PDF report + PowerPoint slides

## Total Tests Passing
~231

## What Is NOT Built Yet
- Phase 2: Quarto report templates ← NEXT
- Phase 2: python-pptx slide deck generator
- Phase 2: auto-methods text + figure captioning
- scVI batch correction → deferred to Phase 6
- ScType-py + SingleR-py annotation
- ADT QC + CLR normalization
- scATAC module (Phase 4)
- Spatial module (Phase 5)
- Multiome module (Phase 6)
- Streamlit UI (Phase 7)
- CLI (Phase 7)

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
- data/processed/GSE194122_cite_harmony.h5ad
- data/processed/GSE194122_cite_harmony_clustered.h5ad
- data/processed/GSE194122_cite_pseudobulk_deg.h5ad

## New Output Structure (run_pipeline.py)
- data/processed/<dataset_id>/01_qc.h5ad → 10_pseudobulk_deg.h5ad
- reports/<dataset_id>/01_qc_report.html → 10_pseudobulk_deg_report.html
