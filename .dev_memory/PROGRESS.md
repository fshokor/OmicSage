# OmicSage — Progress
> Last updated: 2026-05-15 (Phase 3 Session 0 complete)

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

## Phase 3 — AI Layer ← CURRENT PHASE
- [x] Session 0 — Shared infrastructure ✅
  - [x] ai/_base.py — AiResult base dataclass
  - [x] ai/_config_gate.py — three-level check_ai_enabled()
  - [x] ai/_audit_log.py — write_audit_record() → logs/llm/*.jsonl
  - [x] ai/_llm_client.py — BioChatter wrapper, provider routing
  - [x] ai/_skill_loader.py — YAML skill loader (built prior session)
  - [x] ai/skills/cluster_annotator.yaml — reference skill pattern
  - [x] config/study_context_template.yaml
  - [x] tests/test_ai_infrastructure.py — 20 tests passing
  - [x] biochatter==0.14.2 installed and verified
  - [x] test_phase0_structure.py encoding fix (utf-8)
- [x] Session 1 — A1: Pipeline advisor 
  - [x] ai/skills/pipeline_advisor.yaml
  - [x] ai/pipeline_advisor.py
  - [x] tests/test_pipeline_advisor.py
- [ ] Session 2 — A2: Clustering advisor (first PubMed RAG use) ← NEXT
  - [ ] ai/skills/clustering_advisor.yaml
  - [ ] ai/clustering_advisor.py
  - [ ] tests/test_clustering_advisor.py
- [ ] Session 3 — B1: Cluster annotator
  - [ ] ai/skills/cluster_annotator.yaml (already drafted — refine here)
  - [ ] ai/cluster_annotator.py
  - [ ] tests/test_cluster_annotator.py
- [ ] Session 4 — B2: DEG validator + literature linker
  - [ ] ai/skills/deg_validator.yaml
  - [ ] ai/deg_validator.py
  - [ ] tests/test_deg_validator.py
- [ ] Session 5 — B3: Coherence reviewer
  - [ ] ai/skills/coherence_reviewer.yaml
  - [ ] ai/coherence_reviewer.py
  - [ ] tests/test_coherence_reviewer.py
  - [ ] reports/<dataset>/analysis_summary.json (schema designed here)
- [ ] Session 6 — A3: Downstream analysis suggester
  - [ ] ai/skills/downstream_suggester.yaml
  - [ ] ai/downstream_suggester.py
  - [ ] tests/test_downstream_suggester.py
- [ ] Session 7 — C1: Narrative generator
  - [ ] ai/skills/narrative_generator.yaml
  - [ ] ai/narrative_generator.py
  - [ ] tests/test_narrative_generator.py
  - [ ] reports/<dataset>/ai_narrative.md
- [ ] Session 8 — C2: Full report + PowerPoint
  - [ ] ai/skills/report_writer.yaml
  - [ ] ai/report_writer.py
  - [ ] tests/test_report_writer.py
- [ ] Session 9 — Milestone validation
  - [ ] tests/test_groundedness.py
  - [ ] End-to-end GSE166635 with ai_features: true
  - [ ] End-to-end GSE166635 with ai_features: false
  - [ ] MILESTONE: groundedness score ≥ 0.85 ← PHASE 3 COMPLETE

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
