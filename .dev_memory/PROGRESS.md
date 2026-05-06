# OmicSage — Progress
> Last updated: 2026-05-06

## Phase 0 — Foundation
- [x] GitHub repo created
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [x] Full repo structure scaffolded
- [x] Docker base images (Python + R) — defined
- [x] CI/CD via GitHub Actions
- [x] YAML config schema defined
- [x] Benchmark dataset downloaded (GSE194122 CITE-seq + multiome)
- [x] .dev_memory/ system initialized and filled
- [x] README with project description and badges

## Phase 1 — Core scRNA Pipeline
- [x] Data ingestion: 10x MEX, H5, AnnData auto-detection
- [x] QC module: MT%, genes/cell, SoupX ambient RNA, Scrublet doublets
- [x] Modality-aware QC: MuData return, GEX/ADT/ATAC splitting
- [x] Normalization: CP10K + log1p + HVG selection (seurat_v3)
- [x] batch_key support for HVG selection
- [x] layers['counts'] + layers['logcounts'] (Seurat convention)
- [x] Normalization report (HTML, self-contained)
- [x] PCA + UMAP + neighbors (reduce.py)  ← COMPLETED THIS SESSION
- [x] Data-driven PC selection (elbow via kneed + variance fallback)
- [x] Dimensionality reduction report (HTML, self-contained)
- [x] Notebook Step 3 — dimensionality reduction section
- [ ] Leiden clustering with resolution sweep  ← NEXT
- [ ] Silhouette score-based resolution selection
- [ ] Harmony + scVI batch correction
- [ ] SingleR annotation
- [ ] DEG: Wilcoxon + pseudobulk
- [ ] GSEA: GO/KEGG/Reactome
- [ ] MILESTONE: Reproduce key findings of Wang et al. 2025 HCC paper

## Phase 2 — Report Engine
- [ ] Quarto QC report template
- [ ] Quarto analysis report template
- [ ] python-pptx slide deck generator
- [ ] Auto-figure captioning
- [ ] Auto-methods text
- [ ] MILESTONE: Biologist receives complete PDF + slides from one command

## Phase 3 — AI Layer
- [ ] BioChatter integration
- [ ] QC threshold suggestion
- [ ] Cluster interpretation
- [ ] PubMed RAG
- [ ] Narrative generator
- [ ] AI audit log
- [ ] Multi-LLM support
- [ ] MILESTONE: AI narrative groundedness score > 0.85

## Phase 4 — scATAC Module
- [ ] Fragment file ingestion and QC
- [ ] Peak calling (MACS3)
- [ ] LSI + ArchR/Signac clustering
- [ ] Motif enrichment + chromVAR
- [ ] Gene activity scores
- [ ] scATAC report template

## Phase 5 — Spatial Module
- [ ] Visium data ingestion
- [ ] Spatially variable genes
- [ ] BayesSpace clustering
- [ ] RCTD deconvolution
- [ ] Spatial report template

## Phase 6 — Multiome Integration
- [ ] WNN joint embedding
- [ ] MOFA+ integration
- [ ] SCENIC+ GRN inference
- [ ] Joint report template

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
