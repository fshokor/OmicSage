# OmicSage — Current Status
> Last updated: June 2026

## Phase
Phase 8 — Docker Containerization (next session)
Phase 7 — Spatial Transcriptomics — FULLY COMPLETE (Sessions 1–5 + Session A + Session B)

## Test Count
~1476 passing, ~58 skipped (as of end of Phase 7 Session B)
Breakdown: 1429 baseline + 20 (test_spatial_impute) + 27 (test_spatial_ingest_formats)

---

## ✅ Phase 0 — Foundation (COMPLETE)
- Repo structure scaffolded
- Docker base images (Python + R) — defined, not yet built locally
- CI/CD via GitHub Actions — configured
- YAML config schema — defined
- .dev_memory/ system — initialized

---

## ✅ Phase 1 — Core scRNA Pipeline (COMPLETE)

### Ingestion (pipeline/modules/qc/ingest.py)
- Auto-detects 10x MEX, H5, AnnData formats
- Moves normalized values out of .X, puts raw counts in .X
- Handles GSE194122 CITE-seq and multiome, GSE166635 HCC
- load_dataset_dir() — scans parent folder, loads all MTX subfolders,
  concatenates with obs['sample'] + obs['batch'] per subfolder name

### QC (pipeline/modules/qc/qc.py)
- Modality-aware: auto-detects GEX / ADT / ATAC from var['feature_types']
- Returns MuData with mdata["rna"], mdata["adt"], mdata["atac"] slots
- MT%, ribo%, hb%, genes/cell, Scrublet doublets
- 56 tests passing in tests/test_qc.py

### Normalization (pipeline/modules/qc/normalize.py)
- CP10K → log1p, layers['counts'] + layers['logcounts']
- HVG selection — top 2000, seurat_v3 flavor, batch_key support
- 12 tests passing in tests/test_normalize.py

### Dimensionality Reduction (pipeline/modules/qc/reduce.py)
- PCA (HVG subset, n_comps=50, arpack), elbow detection via kneed
- Neighbor graph, UMAP (always), t-SNE (optional)
- 12 tests passing in tests/test_reduce.py

### Leiden Clustering (pipeline/modules/clustering/cluster.py)
- Resolution sweep, silhouette scoring, best_resolution_override
- obs['leiden'] + obs['leiden_harmony'] convenience keys
- compute_ari() for cross-annotation comparison
- 16 tests passing in tests/test_cluster.py

### Cell-Type Annotation (pipeline/modules/annotation/annotate.py)
- CellTypist + ScType + SingleR + scANVI — 4-way weighted vote
- obs['cell_type_vote'], obs['cell_type_confidence']
- 18 tests passing, 1 skipped in tests/test_annotate.py

### DEG (pipeline/modules/downstream/deg.py)
- Wilcoxon rank-sum, one-vs-rest, BH correction
- exclude_gene_prefixes param
- 11 tests passing in tests/test_deg.py

### GSEA (pipeline/modules/downstream/gsea.py)
- ORA via gseapy.enrichr — GO BP / KEGG / Reactome
- direction param: up / down / both
- 8 tests passing in tests/test_gsea.py (Enrichr mocked)

### Harmony Batch Correction (pipeline/modules/integration/harmony_correct.py)
- obsm['X_pca_harmony'], obsm['X_umap_harmony']
- neighbors_harmony graph in uns + obsp
- 13 tests passing in tests/test_harmony.py

### Pseudobulk DEG (pipeline/modules/downstream/pseudobulk_deg.py)
- DESeq2 Wald tests via pydeseq2, one-vs-rest per cell type
- min_cells + min_samples filters with graceful skip
- 14 tests passing in tests/test_pseudobulk_deg.py

### Reports (reports/)
- normalization_report.py, reduce_report.py, cluster_report.py,
  annotate_report.py, deg_report.py, gsea_report.py,
  harmony_report.py, pseudobulk_deg_report.py
- combined_report.py → 00_combined_report.html (tabbed)

### Pipeline Runner (run_pipeline.py)
- Config-driven, --from-step / --to-step / --step / --force flags
- Step-level .h5ad caching, startup validation

### Config System (config/)
- config/schema.yaml — platform master schema
- config/runs/GSE194122.yaml
- config/runs/GSE166635.yaml
- config/runs/GSE194122_multiome.yaml

### MILESTONE: Wang et al. 2025 HCC Benchmark ✓
- Full pipeline on GSE166635, key markers and pathways recovered

---

## ✅ Phase 2 — Report Engine (COMPLETE)
- Per-step HTML reports for all Phase 1 steps
- combined_report.py producing tabbed 00_combined_report.html
- MILESTONE: biologist receives full report from one command ✓

---

## ✅ Phase 3 — AI Layer (COMPLETE / ⏸ PAUSED BY DECISION)
Decision (2026-05-17): AI layer fully built and tested. ai_features: false
is the default going forward. Code is NOT deleted. See DECISIONS.md.

Modules built (all tested):
- ai/pipeline_advisor.py (A1) — 13 tests
- ai/clustering_advisor.py (A2) — 22 tests
- ai/cluster_annotator.py (A3) — 23 tests
- ai/deg_validator.py (B1) — 25 tests
- ai/coherence_reviewer.py (B2) — 22 tests
- ai/downstream_suggester.py (B3) — 20 tests
- ai/narrative_generator.py (C1) — 18 tests
- ai/report_writer.py (C2) — 20 tests
- ai/report_reviewer.py (D1) — 18 tests
- Infrastructure (_base, _config_gate, _audit_log, _llm_client,
  _skill_loader) — 20 tests

Manual review prompts:
- ai/manual_review/CITE_MASTER_PROMPT.md
- ai/manual_review/SPATIAL_MASTER_PROMPT.md (v1.1 — June 2026)

---

## ✅ Phase 4 — CITE-seq Pipeline (COMPLETE)

### Modules
- adt_normalize.py — CLR normalization via muon
- adt_doublets.py — marker-pair scoring
- adt_reduce.py — PCA + UMAP on CLR matrix (obsm: X_pca_adt, X_umap_adt)
- adt_harmony.py — harmonypy direct (runtime shape detection)
- adt_annotate.py — Leiden on ADT embedding, obs['adt_celltype']
- cite_integration.py — MOFA+ + totalVI (WNN deferred)
- cite_deg.py — Wilcoxon DPE + cross-modal RNA DEG
- cite_gsea.py — ORA on RNA arm with ADT labels
- cite_corr.py — within-cell-type Spearman r + Fisher z-transform
- cite_epitope.py — protein-level cell type marker summary

### Runner + Reports
- run_cite_pipeline.py (cite_01–cite_10 steps)
- All cite reports + cite_combined_report.py → cite_00_combined_report.html

### Benchmark Results (GSE194122 BMMC CITE-seq, 21,778 cells)
- totalVI recommended embedding: iLISI=0.229, cLISI=1.000, ASW=0.712

---

## ✅ Phase 5 — Multiome RNA + ATAC (COMPLETE)

### Modules
- atac_qc.py — TSS enrichment, nucleosome signal, fragment-based QC
- atac_reduce.py — TF-IDF → LSI → UMAP (obsm: X_lsi, X_umap_atac)
- atac_annotate.py — gene activity scores + RNA-guided annotation
- multiome_integration.py — MultiVI joint embedding (obsm: X_multivi)
- multiome_deg.py — joint RNA DEG + differential accessibility
- run_multiome_pipeline.py + combined report

---

## ✅ Phase 6 — GRN Inference (COMPLETE)

### Modules
- grn_infer.py — pyscenic + decoupler stack
- grn_report.py
- run_grn_pipeline.py

---

## ✅ Phase 7 — Spatial Transcriptomics (FULLY COMPLETE)
> Sessions 1–5 (core pipeline) + Session A (imputation) + Session B (new formats)

### Pipeline Modules (pipeline/modules/spatial/)

#### spatial_ingest.py
- Unified dispatcher with auto-detection
- Supported: visium | h5ad | benchmark | visium_hd | xenium
- visium_hd: via spatialdata_io.visium_hd(), bin_size param (2/8/16µm)
- xenium: via spatialdata_io.xenium(), cell-level, obsm["spatial"] from centroids
- ENSEMBL ID support, MT gene handling, RGBA alpha-stripping
- library_key auto-detection from ingest provenance

#### spatial_qc.py
- QC metrics: total_counts, n_genes_by_counts, pct_counts_mt
- Technology-aware thresholds documented (Visium/Visium HD/Xenium)
- filter_spots param, summary stats in provenance

#### spatial_reduce.py
- normalize_total → log1p → HVG → PCA → spatial neighbours graph
- flavor="cell_ranger" for ATAC; flavor="seurat" for RNA
- sq.gr.spatial_neighbors with coord_type auto-detect

#### spatial_cluster.py
- Leiden on transcriptomic KNN (NOT spatial graph)
- Moran's I SVG detection on spatial graph
- uns["moranI"] DataFrame with FDR-corrected p-values

#### spatial_deconvolve.py
- Two-tier: NNLS (default, ~200 MB RAM) + cell2location (opt-in)
- per_sample loop support, identical output contract
- method="none" for Xenium (cell-level, deconvolution not applicable)
- obsm["q05_cell_abundance_w_sf"] canonical key

#### spatial_downstream.py
- Eight analyses: cell-type specific gene expression correlation,
  cell-type specific SVGs, spatial co-occurrence, neighbourhood
  enrichment, ligand-receptor communication (OmniPath), pathway
  enrichment (gseapy prerank on Moran's I), region clustering
- _sanitize_gsea_df() fixes gseapy dtype bug (object → numeric)

#### spatial_impute.py (Session A — NEW)
- Tangram (default, clusters mode) + gimVI (opt-in)
- tangram_mode="clusters": maps cell-type signatures, negligible RAM
- tangram_mode="cells": maps individual cells, subsampled via
  max_cells_per_type (default 500) to prevent OOM
- project_genes fix: cluster_label passed, return value captured
- Storage: float32 numpy array in obsm["imputed_expression"],
  gene names in uns["omicsage_spatial_impute"]["outputs"]["genes_imputed"]
- Graceful skip when sc_reference_path absent or enabled=false

### Reports (reports/templates/spatial/)
- spatial_qc_report.py, spatial_reduce_report.py, spatial_cluster_report.py,
  spatial_deconvolve_report.py, spatial_downstream_report.py
- spatial_impute_report.py (Session A — NEW)
  - Section 2: clusters mode renders explanatory note (not skip)
  - Section 4: log-normalises measured counts before Spearman r
  - library_key passed to sq.pl.spatial_scatter for multi-sample data
- spatial_combined_report.py — 6 tabs (QC/Reduce/Cluster/Deconvolve/Downstream/Impute)
  Lightbox click-to-expand, max-width CSS, panoramic stitch cap (1800px)

### Runner (run_spatial_pipeline.py)
- 7 steps: ingest → qc → reduce → cluster → deconvolve → downstream → impute
- Checkpoints: 01_ingested.h5ad through 07_imputed.h5ad
- --step / --from-step / --to-step / --force flags

### Tests
- tests/test_spatial_ingest_qc.py — 73 tests (updated: visium_hd/xenium
  stub tests replaced with implemented assertions)
- tests/test_spatial_ingest_formats.py — 27 tests (Session B — NEW)
  Visium HD + Xenium loaders, all mocked (no spatialdata-io required in CI)
- tests/test_spatial_impute.py — 20 tests (Session A — NEW)
  Tangram fully mocked via sys.modules patch

### Configs (config/runs/)
- kuppe_heart.yaml — impute block added (tangram, clusters mode)
- visium_hd_mouse_brain.yaml — NEW (10x Mouse Brain, 8µm bins)
- xenium_breast.yaml — NEW (10x Human Breast Cancer, targeted panel)

### Benchmark Data (scripts/download_spatial_benchmark.py — NEW)
- Kuppe Visium + snRNA-seq: auto-download from Figshare
- Visium HD Mouse Brain: manual (10x website, free)
- Xenium Breast Cancer: manual (10x website, account required)

### Documentation
- SPATIAL_MODULE_DOCS.md — fully updated (Sessions 1–5 + A + B)
- ai/manual_review/SPATIAL_MASTER_PROMPT.md — v1.1
  (Task 9 imputation review, Xenium/Visium HD guidance throughout)

### Benchmark Results (Kuppe et al. 2022 Human Heart, 11,706 spots)
- 4 control samples (control_P1/P7/P8/P17)
- 8 spatial clusters identified
- Top SVGs: cardiac sarcomere genes (MYH7, MYH6, TNNT2)
- Deconvolution (NNLS): cardiomyocytes dominant in myocardial zones
- Imputation (Tangram clusters mode): 2,000 genes imputed from 41,663-cell
  snRNA-seq reference

### Key Bugs Fixed in Phase 7
- Python 3.11 f-string backslash syntax errors
- gseapy dtype bug (numeric columns → object) — _sanitize_gsea_df()
- squidpy API: show=False, img_res_key vs img_key, fig_params
- harmonypy ≥2.0 Z_corr shape detection
- Tangram OOM: cells mode → clusters mode default
- Tangram project_genes: cluster_label kwarg + return value capture
- imputed_expression: DataFrame → numpy array (AnnData obsm rejects str/DataFrame)
- pd.read_json pandas ≥2.0: StringIO wrapper required
- spatial_scatter multi-sample: library_key from ingest provenance

---

## 🔄 Phase 8 — Docker Containerization (NEXT SESSION)

### Goal
Working `docker-compose up` that runs the full spatial pipeline on the
Kuppe benchmark inside a container, with results on the host filesystem.

### Deliverables
- docker/Dockerfile.pipeline — conda-based Python analysis environment
- docker/Dockerfile.dev — development image with hot-reload
- docker-compose.yml — volume mounts for data/, reports/, config/
- .dockerignore — exclude data/, checkpoints/, __pycache__
- scripts/run_docker.sh — convenience wrapper

### Key Decisions Already Made
- Base image: continuumio/miniconda3 (matches WSL2 conda env)
- Two-stage build: deps layer (environment.yml) + code layer
- Dev workflow: mount repo as volume, no rebuild on code change
- Entrypoint: python run_spatial_pipeline.py "$@"

---

## ⏳ Phase 9 — Streamlit No-Code UI (PLANNED)

Drag-drop upload, guided config, progress tracking, report viewer.
Calls existing runner and modules directly.

---

## ⏳ Phase 10 — Nextflow DSL2 (PLANNED)

Wraps Python modules as Nextflow processes.
Each step = one process block. Python modules unchanged.

---

## What Is NOT Built Yet
- Streamlit UI (Phase 9)
- Nextflow DSL2 wrapping (Phase 10)
- CLI (click-based omicsage.py) — deferred
- WNN joint embedding — deferred (pynndescent hang unresolved)
- DSB normalization — deferred (no empty droplets in GSE194122)
- MERFISH / CODEX ingest — stubs only (NotImplementedError)
- scATAC standalone pipeline — absorbed into Multiome phase
- Docker image — built next session

---

## Dependencies Added in Phase 7 (add to environment.yml + requirements-ci.txt)
- tangram-sc (pip only, no conda package)
- spatialdata-io (pip: pip install spatialdata-io)

---

## Processed Data Files

### Spatial (Kuppe et al. 2022 Human Heart)
- data/processed/kuppe_heart/01_ingested.h5ad
- data/processed/kuppe_heart/02_qc.h5ad
- data/processed/kuppe_heart/03_reduced.h5ad
- data/processed/kuppe_heart/04_clustered.h5ad
- data/processed/kuppe_heart/05_deconvolved.h5ad
- data/processed/kuppe_heart/06_downstream.h5ad
- data/processed/kuppe_heart/07_imputed.h5ad
- reports/kuppe_heart/00_spatial_combined_report.html

### CITE-seq (GSE194122)
- data/processed/GSE194122/cite_01 through cite_10 .h5mu files

### scRNA (GSE166635 HCC benchmark)
- data/processed/GSE166635/05_annotated.h5ad
