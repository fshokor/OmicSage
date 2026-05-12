# OmicSage — Progress
> Last updated: 2026-05-12 (session 9 — Phase 2 complete, Phase 3 starting)

## Phase 0 — Foundation
- [x] GitHub repo created
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [x] Full repo structure scaffolded
- [x] Docker base images (Python + R) — defined
- [x] CI/CD via GitHub Actions
- [x] YAML config schema defined (config/schema.yaml — updated this session)
- [x] Benchmark dataset downloaded (GSE194122 CITE-seq + multiome)
- [x] .dev_memory/ system initialized and filled
- [x] README with project description and badges

## Phase 1 — Core scRNA Pipeline
- [x] Data ingestion: 10x MEX, H5, AnnData auto-detection
- [x] Data ingestion: multi-sample MTX parent folder (load_dataset_dir) ← NEW
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
- [x] Generic pipeline runner (run_pipeline.py) ← NEW
- [x] Config system (config/schema.yaml + config/runs/*.yaml) ← NEW
- [x] Per-dataset run configs: GSE194122, GSE166635, GSE194122_multiome ← NEW
- [x] MILESTONE: Reproduce key findings of Wang et al. 2025 HCC paper ← DONE

## Phase 2 — Report Engine ✅ COMPLETE
- [x] Architectural decision: keep HTML, defer Quarto + pptx to Phase 3 AI layer
- [x] Combined tabbed report (reports/combined_report.py)
- [x] Wired into run_pipeline.py — auto-generates at end of every run
- [x] MILESTONE: 00_combined_report.html generated automatically — 7 tabs, tested on GSE166635

## Phase 3 — AI Layer ← CURRENT PHASE
- [ ] BioChatter integration ← NEXT
- [ ] QC threshold suggester (ai/threshold_suggester.py)
- [ ] Cluster interpreter (ai/cluster_interpreter.py)
- [ ] PubMed RAG tied to DEG results
- [ ] Narrative generator for combined report (+ PDF/slides output)
- [ ] AI audit log (logs/llm/ JSONL)
- [ ] Multi-LLM support (Claude / GPT-4o / local Ollama)
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
- [ ] scVI / MultiVI batch correction
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
