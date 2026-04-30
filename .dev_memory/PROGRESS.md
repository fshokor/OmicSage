# OmicSage — Progress Checklist
Tick items as they are completed. Do not remove items — mark [x] done.

---

## Phase 0 — Foundation (Week 1-2)

- [x] GitHub repo created (https://github.com/fshokor/OmicSage)
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [x] Full repo structure scaffolded (all dirs and stubs)
- [x] Docker base image — Python (Dockerfile.python defined)
- [x] Docker base image — R (Dockerfile.r defined)
- [x] Docker base images BUILT locally ← **TODO: run docker build**
- [x] CI/CD via GitHub Actions (ci.yml + docker-publish.yml)
- [x] YAML config schema defined (config/schema.yaml)
- [ ] Benchmark dataset downloaded (GSE166635 HCC scRNA-seq)
- [x] .dev_memory/ system initialized and all 5 files filled
- [ ] README.md with project description and badges ← **TODO next session**
- [ ] .gitignore created ← **TODO next session**
- [ ] DVC initialized (`dvc init`)
- [ ] First commit pushed to GitHub
- [ ] CI passes on GitHub (green checkmark)

---

## Phase 1 — Core scRNA Pipeline (Week 2-6)

- [ ] Data ingestion: 10x MEX, H5, AnnData auto-detection
- [ ] QC module: MT%, genes/cell, SoupX ambient RNA, Scrublet doublets
- [ ] Normalization + HVG selection (scran)
- [ ] PCA + UMAP + t-SNE
- [ ] Harmony + scVI batch correction
- [ ] Leiden clustering with resolution sweep
- [ ] SingleR cell type annotation
- [ ] DEG: Wilcoxon
- [ ] DEG: pseudobulk DESeq2
- [ ] GSEA: GO/KEGG/Reactome
- [ ] **MILESTONE**: Reproduce Wang et al. 2025 HCC key findings from raw counts

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
