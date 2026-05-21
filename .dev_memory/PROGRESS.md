# OmicSage — Progress
> Last updated: 2026-05-17

## Phase 0 — Foundation ✅ COMPLETE
- [x] GitHub repo created
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [x] Full repo structure scaffolded
- [x] Docker base images (Python + R) — defined
- [x] CI/CD via GitHub Actions
- [x] YAML config schema defined (config/schema.yaml)
- [x] Benchmark dataset downloaded (GSE194122 CITE-seq + multiome)
- [x] .dev_memory/ system initialized and filled
- [x] README with project description and badges

## Phase 1 — Core scRNA Pipeline ✅ COMPLETE
- [x] Data ingestion: 10x MEX, H5, AnnData auto-detection
- [x] Data ingestion: multi-sample MTX parent folder (load_dataset_dir)
- [x] QC module: MT%, genes/cell, Scrublet doublets
- [x] Modality-aware QC: MuData return, GEX/ADT/ATAC splitting
- [x] Normalization: CP10K + log1p + HVG selection (seurat_v3)
- [x] batch_key support for HVG selection
- [x] layers['counts'] + layers['logcounts'] (Seurat convention)
- [x] Normalization report (HTML, self-contained)
- [x] PCA + UMAP + neighbors (reduce.py)
- [x] Data-driven PC selection (elbow via kneed + variance fallback)
- [x] Dimensionality reduction report (HTML, self-contained)
- [x] Leiden clustering with resolution sweep (cluster.py)
- [x] Silhouette score per resolution + best_resolution_override
- [x] Clustering report (HTML, self-contained)
- [x] Cell-type annotation: CellTypist + marker scoring + majority vote
- [x] Annotation report (HTML, self-contained)
- [x] DEG: Wilcoxon rank-sum, one-vs-rest, BH correction
- [x] DEG report: volcano plots + dot plot + summary table
- [x] GSEA: ORA via gseapy.enrichr — GO BP / KEGG / Reactome
- [x] GSEA report: bar charts + bubble plot + direction badges
- [x] Harmony batch correction
- [x] Harmony report: batch composition, mixing metrics, PC shift, UMAP comparison
- [x] Clustering on Harmony-corrected embedding
- [x] Pseudobulk DEG: pydeseq2-based, one-vs-rest per cell type
- [x] Pseudobulk DEG report: skipped groups section
- [x] Notebook Steps 1–10 complete
- [x] Generic pipeline runner (run_pipeline.py)
- [x] Config system (config/schema.yaml + config/runs/*.yaml)
- [x] Per-dataset run configs: GSE194122, GSE166635, GSE194122_multiome
- [x] MILESTONE: Reproduce key findings of Wang et al. 2025 HCC paper ✅

## Phase 2 — Report Engine ✅ COMPLETE
- [x] Architectural decision: keep HTML, defer Quarto + pptx to Phase 3 AI layer
- [x] Combined tabbed report (reports/combined_report.py)
- [x] Wired into run_pipeline.py — auto-generates at end of every run
- [x] MILESTONE: 00_combined_report.html generated automatically — 7 tabs ✅

## Phase 3 — AI Layer ✅ BUILT / ⏸ PAUSED BY DECISION (2026-05-17)
> All modules built and tested (466+ tests passing). Development stopped.
> Manual pipeline is the primary path. ai_features: false is the default.
> See DECISIONS.md for full rationale. Code is NOT deleted — stays intact.

- [x] Session 0 — Shared infrastructure (20 tests)
- [x] Session 1 — A1: Pipeline advisor (13 tests)
- [x] Session 2 — A2: Clustering advisor (22 tests)
- [x] Session 3 — B1: Cluster annotator (23 tests)
- [x] Session 4 — B2: DEG validator + literature linker (25 tests)
- [x] Session 5 — B3: Coherence reviewer (22 tests)
- [x] Session 6 — A3: Downstream analysis suggester (20 tests)
- [x] Session 7 — C1: Narrative generator (18 tests)
- [x] Session 8 — C2: Full report + PowerPoint (20 tests)
- [x] Session 9 — D1: Report reviewer (18 tests)
- [x] Session 10 — Milestone validation complete

## Phase 1 — Annotation Module 
- [x] Annotated h5ad output with cell_type column (SingleR-based)
- [x] Annotation section in 00_combined_report.html
- [x] MILESTONE: Full Phase 1 pipeline with annotation on GSE166635

## Phase 4 — CITE-seq Module 🔄
- [ ] ADT normalisation (CLR + optional DSB) — `pipeline/modules/cite/adt_normalize.py`
- [ ] ADT dimensionality reduction (PCA on CLR-normalised ADT)
- [ ] WNN joint embedding (RNA + ADT)
- [ ] WNN UMAP + cluster report
- [ ] Protein-level annotation validation (CD303, CD138, CD34, CD14/CD16)
- [ ] CITE-seq report tab
- [ ] **MILESTONE**: WNN UMAP showing RNA + protein integrated embedding

## Phase 5 — Multiome Integration (RNA + ATAC)
- [ ] scVI / MultiVI batch correction
- [ ] RNA + ATAC joint ingestion and QC
- [ ] LSI dimensionality reduction for ATAC
- [ ] WNN joint embedding (RNA + ATAC)
- [ ] MOFA+ integration
- [ ] Peak calling (MACS3)
- [ ] Motif enrichment (chromVAR)
- [ ] Gene activity scores
- [ ] SCENIC+ GRN inference
- [ ] Joint report template
- NOTE: scATAC-seq standalone is NOT a separate phase —
        ATAC is always analysed jointly with RNA in Multiome

## Phase 6 — Spatial Module
- [ ] Visium data ingestion (10x Space Ranger output)
- [ ] Spatially variable genes (squidpy)
- [ ] BayesSpace clustering
- [ ] RCTD deconvolution
- [ ] Spatial report template
- [ ] Support for MERFISH and Xenium (future)

## Phase 7 — User Interfaces
- [ ] Streamlit web UI
- [ ] Click CLI
- [ ] Project dashboard
- [ ] Shared protocol library

## Phase 8 — Benchmarking + Paper
- [ ] Benchmark on 5 published datasets
- [ ] User study
- [ ] Paper submitted
- [ ] bioRxiv preprint
- [ ] v1.0 public release
- [ ] omicsage.io website launched
