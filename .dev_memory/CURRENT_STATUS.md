# OmicSage — Current Status
> Updated: 2026-06-17
> Phase: Nextflow DSL2 complete → Tutorial / Demo / Colab next

---

## What Is Built and Tested

### Phase 1 — scRNA Pipeline ✅
Full pipeline: ingest → QC → normalize → reduce → cluster → annotate →
DEG → GSEA → Harmony → cluster_harmony → pseudobulk DEG.
Config-driven runner. Combined tabbed HTML report.

### Phase 2 — Report Engine ✅
Combined tabbed HTML report (`00_combined_report.html`).
Auto-generated after every pipeline run.

### Phase 3 — AI Layer ✅ (built, paused by decision)
BioChatter-based. 8 modules: Pipeline Advisor, Clustering Advisor,
Cluster Annotator, DEG Validator, Coherence Reviewer, Downstream Suggester,
Narrative Generator, Report Writer + PPTX.
ai_features: false is the default. Code intact, not deleted.

### Phase 4 — CITE-seq Pipeline ✅
Full pipeline: ADT normalize → doublets → reduce → harmony → annotate →
integration → DEG/DPE → GSEA → correlation → epitope characterization.
Combined report + runner (run_cite_pipeline.py).

### Phase 5 — Multiome Pipeline ✅
ATAC QC → LSI reduce → annotate → MultiVI integration →
multiome DEG → GRN with pyscenic + decoupler.
Runner + combined report (run_multiome_pipeline.py).

### Phase 6 — GRN Inference ✅
GRN inference module complete and integrated into multiome pipeline.

### Phase 7 — Spatial Transcriptomics Pipeline ✅
Full pipeline: ingest → QC → reduce → cluster → deconvolve (NNLS/cell2location)
→ downstream → impute (Tangram/gimVI).
Configs: kuppe_heart.yaml, visium_hd_mouse_brain.yaml, xenium_breast.yaml.
Runner + combined report (run_spatial_pipeline.py).

### Phase 8 — Docker Containerization ✅
Single image (omicsage:latest, ~7.6GB) covers all four modalities.
Entrypoint routes by OMICSAGE_MODALITY env var or first CLI argument.
docker-compose.yml with 5 services. scripts/run_docker.sh convenience wrapper.

### Phase 9 — Streamlit UI ✅
4-page UI: Dataset → Configure → Run → Report.
Python runner only (no Docker overhead in UI).
GPU toggle via OMICSAGE_GPU=1 env var.
Live log streaming via background thread + queue.
Save/load config YAML. Project history sidebar.

### Phase 10 — Nextflow DSL2 Orchestration ✅ ← COMPLETED THIS SESSION
Full DSL2 layer for all 4 modalities.

**Files delivered:**
```
main.nf                              ← entry point, routes by --modality
nextflow.config                      ← Docker mount, GPU flag, memory, profiles
pipeline/workflows/
  scrna.nf                           ← 10-step scRNA workflow
  cite.nf                            ← 10-step CITE-seq workflow
  multiome.nf                        ← 6-step Multiome workflow
  spatial.nf                         ← 7-step Spatial workflow
pipeline/modules/
  scrna/   (10 modules)              ← qc → pseudobulk
  cite/    (10 modules)              ← normalize_adt → epitope_characterisation
  multiome/(6 modules)               ← atac_qc → multiome_grn
  spatial/ (7 modules)               ← ingest → impute
```

**Key design decisions:**
- Each module calls Python runner with full Python path (entrypoint bypassed via --entrypoint '')
- Config passed as val (not path) to avoid Nextflow staging stripping the path
- Skip guard in every module: reads YAML, exits 0 + touches sentinel if enabled: false
- Streamlit UI uses Python runner only (step control, --from-step/--to-step)
- Nextflow is CLI-only for HPC/cloud full pipeline runs
- GPU: --gpu true flag in nextflow.config passes --gpus all + OMICSAGE_GPU=1

**How to run via Nextflow (CLI/HPC):**
```bash
# scRNA full pipeline
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile local \
  -ansi-log false

# Resume after crash
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile local \
  -resume

# With GPU
nextflow run main.nf \
  --config config/runs/GSE194122_cite.yaml \
  --modality cite \
  -profile local \
  --gpu true

# HPC (SLURM)
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile slurm
```

---

## Current Test Count
**1472 passing, 58 skipped** (as of 2026-06-16)

---

## What Is NOT Yet Built
- Tutorial notebook (GSE166635 walkthrough) ← NEXT
- Demo video (Streamlit screen recording)
- Google Colab notebook (GPU-accessible, Drive-mounted)
- README rewrite (what's built + how to use it)
- CLI (Click-based project manager)
- Paper / benchmarking
- omicsage.io website
