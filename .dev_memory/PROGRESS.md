# OmicSage — Progress Checklist
Tick items as they are completed. Do not remove items — mark [x] done.

---

## Phase 0 — Foundation ✅ COMPLETE
- [x] GitHub repo created (https://github.com/fshokor/OmicSage)
- [x] Name: OmicSage
- [x] Domain confirmed: omicsage.io
- [x] Full repo structure scaffolded
- [x] Dockerfile.python defined + DEBIAN_FRONTEND fix applied
- [x] Dockerfile.r defined + DEBIAN_FRONTEND fix applied
- [x] CI/CD: ci.yml (lint + pytest + docker build main-only)
- [x] CI/CD: docker-publish.yml (GHCR on version tags)
- [x] config/schema.yaml complete
- [x] .dev_memory/ all 5 files initialized
- [x] README.md, LICENSE, .gitignore
- [x] 52/52 tests passing locally
- [x] Lint job green on GitHub Actions
- [x] Python Tests job green on GitHub Actions
- [x] Docker build fix: DEBIAN_FRONTEND=noninteractive in both images
- [x] .gitkeep files added to all empty dirs
- [ ] Docker images built locally (run manually — too heavy for CI)
- [ ] GSE166635 benchmark data downloaded
- [ ] DVC initialized

## Phase 1 — Core scRNA Pipeline (Week 2-6) ← WE ARE HERE
- [ ] Data ingestion module (10x MEX, H5, AnnData auto-detection)
- [ ] QC module: MT%, genes/cell, Scrublet doublets, SoupX ambient RNA
- [ ] Normalization + HVG selection (scran)
- [ ] PCA + UMAP + t-SNE
- [ ] Harmony + scVI batch correction
- [ ] Leiden clustering with resolution sweep
- [ ] SingleR cell type annotation
- [ ] DEG: Wilcoxon
- [ ] DEG: pseudobulk DESeq2
- [ ] GSEA: GO/KEGG/Reactome
- [ ] MILESTONE: Reproduce Wang et al. 2025 HCC key findings from raw counts

---

## Phase 2 — Report Engine (Week 6-9)

- [ ] Quarto QC report template
- [ ] Quarto analysis report template (clustering, DEG, pathway)
- [ ] python-pptx slide deck generator
- [ ] Auto-figure captioning from metadata
- [ ] Auto-methods text from config + software versions
- [ ] **MILESTONE**: Complete PDF + slides from one command, no coding

---

## Phase 3 — AI Layer (Week 9-13)

- [ ] BioChatter integration (ai/biochatter_client.py implemented)
- [ ] QC threshold suggestion via augmented prompts
- [ ] Cluster interpretation (markers → LLM → cell type + evidence)
- [ ] PubMed RAG tied to DEG results
- [ ] Biological narrative generator for reports
- [ ] AI audit log (all LLM calls logged)
- [ ] Multi-LLM support (Claude, GPT-4o, local Ollama)
- [ ] **MILESTONE**: AI report narrative groundedness score > 0.85

---

## Phase 4 — scATAC Module (Week 13-16)

- [ ] Fragment file ingestion and QC
- [ ] Peak calling (MACS3)
- [ ] LSI + ArchR/Signac clustering
- [ ] Motif enrichment + chromVAR
- [ ] Gene activity scores
- [ ] scATAC report template

---

## Phase 5 — Spatial Module (Week 16-19)

- [ ] Visium data ingestion
- [ ] Spatially variable genes
- [ ] BayesSpace clustering
- [ ] RCTD deconvolution
- [ ] Spatial report template

---

## Phase 6 — Multiome Integration (Week 19-22)

- [ ] WNN joint embedding
- [ ] MOFA+ integration
- [ ] SCENIC+ GRN inference
- [ ] Joint report template

---

## Phase 7 — User Interfaces (Week 22-25)

- [ ] Streamlit web UI (drag-drop upload, guided config, progress)
- [ ] Click CLI polished (create-project, run, list, compare)
- [ ] Project dashboard
- [ ] Shared protocol library

---

## Phase 8 — Benchmarking + Paper (Week 25-30)

- [ ] Benchmark on 5 published datasets
- [ ] User study (biologists + bioinformaticians)
- [ ] Paper written (target: Nature Methods or Bioinformatics)
- [ ] bioRxiv preprint posted
- [ ] v1.0 public release
- [ ] omicsage.io website launched
