# OmicSage — Progress
> Last updated: 2026-06-17

## Phase 0 — Foundation ✅ COMPLETE
- [x] GitHub repo created
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [x] Full repo structure scaffolded
- [x] CI/CD via GitHub Actions
- [x] YAML config schema defined (config/schema.yaml)
- [x] Benchmark datasets downloaded (GSE194122, GSE166635, kuppe_heart)
- [x] .dev_memory/ system initialized and filled
- [x] README with project description and badges

## Phase 1 — Core scRNA Pipeline ✅ COMPLETE
- [x] Data ingestion: 10x MEX, H5, AnnData, multi-sample MTX
- [x] QC: MT%, genes/cell, Scrublet doublets
- [x] Normalization: CP10K + log1p + HVG selection
- [x] PCA + UMAP + neighbors
- [x] Leiden clustering with resolution sweep + silhouette scoring
- [x] Cell-type annotation: CellTypist + marker scoring + vote
- [x] DEG: Wilcoxon rank-sum, BH correction
- [x] GSEA: ORA via gseapy — GO BP / KEGG / Reactome
- [x] Harmony batch correction + re-clustering
- [x] Pseudobulk DEG (pydeseq2)
- [x] Generic pipeline runner (run_scrna_pipeline.py)
- [x] Config system (config/schema.yaml + config/runs/*.yaml)
- [x] MILESTONE: Reproduce key findings of Wang et al. 2025 HCC paper ✅

## Phase 2 — Report Engine ✅ COMPLETE
- [x] Combined tabbed HTML report (reports/combined_report.py)
- [x] Wired into run_pipeline.py — auto-generates at end of every run
- [x] MILESTONE: 00_combined_report.html — 7 tabs, tested on GSE166635 ✅

## Phase 3 — AI Layer ✅ BUILT / ⏸ PAUSED BY DECISION (2026-05-17)
- [x] 8 modules built and tested (466+ tests)
- [x] ai_features: false default — pipeline works without AI
- [x] Code intact, not deleted

## Phase 4 — CITE-seq Pipeline ✅ COMPLETE
- [x] ADT normalization, doublets, reduce, harmony, annotate
- [x] RNA + ADT integration (WNN / totalVI / MOFA+)
- [x] DEG/DPE, GSEA, protein-RNA correlation, epitope characterization
- [x] Combined report + runner (run_cite_pipeline.py)

## Phase 5 — Multiome Pipeline ✅ COMPLETE
- [x] ATAC QC, LSI reduce, annotate (label transfer)
- [x] MultiVI integration
- [x] Multiome DEG
- [x] GRN inference (pyscenic + decoupler)
- [x] Runner + combined report (run_multiome_pipeline.py)

## Phase 6 — GRN Inference ✅ COMPLETE
- [x] GRN module integrated into multiome pipeline

## Phase 7 — Spatial Transcriptomics Pipeline ✅ COMPLETE
- [x] Ingest: Visium, H5AD, Visium HD, Xenium (spatialdata-io)
- [x] QC, reduce, cluster
- [x] Deconvolution: NNLS (default) + cell2location (opt-in)
- [x] Downstream: SVG, niche analysis, CCC (squidpy)
- [x] Imputation: Tangram + gimVI
- [x] Combined report + runner (run_spatial_pipeline.py)
- [x] Configs: kuppe_heart, visium_hd_mouse_brain, xenium_breast

## Phase 8 — Docker Containerization ✅ COMPLETE
- [x] Single image (omicsage:latest) all 4 modalities
- [x] Entrypoint routing by modality
- [x] docker-compose.yml (5 services)
- [x] scripts/run_docker.sh convenience wrapper
- [x] GPU opt-in via OMICSAGE_GPU=1

## Phase 9 — Streamlit UI ✅ COMPLETE
- [x] Page 1: Dataset (modality, name, ID, organism, paths)
- [x] Page 2: Configure (step selector, per-step params, all 4 modalities)
- [x] Page 3: Run (step range, --force, GPU toggle, live log streaming)
- [x] Page 4: Report (embedded HTML, download button)
- [x] Save/load config YAML
- [x] Project history sidebar
- [x] Python runner only (no Docker overhead in UI)

## Phase 10 — Nextflow DSL2 Orchestration ✅ COMPLETE ← DONE THIS SESSION
- [x] main.nf — routes by --modality
- [x] nextflow.config — Docker mount, GPU, memory, profiles (local/slurm/ci/test)
- [x] 4 workflow files (scrna, cite, multiome, spatial)
- [x] 33 module files with skip guard (enabled: false respected)
- [x] Entrypoint bypass (--entrypoint '')
- [x] GPU flag (--gpu true → --gpus all + OMICSAGE_GPU=1)
- [x] Tested: scRNA pipeline runs end-to-end via Nextflow
- [x] MILESTONE: nextflow run main.nf works for all 4 modalities ✅

## Phase 11 — Portfolio & Outreach ← NEXT
- [ ] Tutorial notebook (GSE166635 scRNA walkthrough, end-to-end)
- [ ] Demo video (Streamlit screen recording, 3-5 min)
- [ ] Google Colab notebook (GPU-accessible, Google Drive mounted)
- [ ] README rewrite (what's built + how to use it, no roadmap)
- [ ] omicsage.io website

## Phase 12 — Benchmarking + Paper
- [ ] Benchmark on 5 published datasets
- [ ] User study (biologists + bioinformaticians)
- [ ] Paper submitted (target: Nature Methods or Bioinformatics)
- [ ] bioRxiv preprint
- [ ] v1.0 public release
