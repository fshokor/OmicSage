# OmicSage — Current Status
> Updated: 2026-06-16
> Phase: Docker containerization complete → Streamlit UI next

---

## What Is Built and Tested

### Phase 1 — scRNA Pipeline ✅
Full pipeline: ingest → QC → normalize → reduce → cluster → annotate →
DEG → GSEA → Harmony → pseudobulk DEG. Config-driven runner.
Combined tabbed HTML report.

### Phase 2 — Report Engine ✅
Combined tabbed HTML report (`00_combined_report.html`).
Auto-generated after every pipeline run.

### Phase 3 — AI Layer ✅
BioChatter-based. 8 modules:
Pipeline Advisor, Clustering Advisor, Cluster Annotator, DEG Validator,
Coherence Reviewer, Downstream Suggester, Narrative Generator,
Report Writer + PPTX. Manual Review Mode (`MASTER_PROMPT.md`).

### Phase 4 — CITE-seq Pipeline ✅
Full pipeline: ADT normalize → doublets → reduce → harmony → annotate →
integration → DEG/DPE → GSEA → correlation → epitope characterization.
Combined report + runner.

### Phase 5 — Multiome Pipeline ✅
ATAC QC → LSI reduce → annotate → MultiVI integration →
GRN with pyscenic + decoupler. Runner + combined report.

### Phase 6 — GRN Inference ✅
GRN inference module complete.

### Phase 7 — Spatial Transcriptomics Pipeline ✅
Full pipeline across Sessions 1–5 + A + B:
- Ingest: Visium / H5AD / benchmark / Visium HD / Xenium (spatialdata-io)
- QC, reduce, cluster, deconvolution (cell2location + NNLS)
- Downstream analyses
- Spatial imputation: Tangram (clusters mode) + gimVI (opt-in)
- Combined report + runner
- Configs: kuppe_heart.yaml, visium_hd_mouse_brain.yaml, xenium_breast.yaml

### Phase 8 — Docker Containerization ✅ ← COMPLETED THIS SESSION
Single image covers all four modalities.

**Files delivered:**
```
docker/Dockerfile.pipeline   ← conda-based, continuumio/miniconda3:23.10.0-1
docker/Dockerfile.dev        ← extends pipeline + Jupyter Lab + pytest tools
docker/entrypoint.sh         ← routes modality via OMICSAGE_MODALITY env var
docker-compose.yml           ← 5 services: scrna, cite, multiome, spatial, dev
.dockerignore                ← excludes data/, __pycache__, .git, reports/output
scripts/run_docker.sh        ← convenience wrapper, handles volumes + GPU flag
```

**How to run any modality:**
```bash
# Via convenience script
./scripts/run_docker.sh scrna --config config/runs/gse166635_hcc.yaml
./scripts/run_docker.sh cite --config config/runs/gse194122_cite.yaml
./scripts/run_docker.sh multiome --config config/runs/gse194122_multiome.yaml
./scripts/run_docker.sh spatial --config config/runs/kuppe_heart.yaml

# Via docker-compose
docker-compose run scrna --config config/runs/gse166635_hcc.yaml
docker-compose run spatial --config config/runs/kuppe_heart.yaml --step ingest

# Dev environment
./scripts/run_docker.sh dev jupyter    # Jupyter Lab on :8888
./scripts/run_docker.sh dev bash       # interactive shell
./scripts/run_docker.sh dev pytest tests/ -q
```

**Key design decisions:**
- One image, four modalities — no image sprawl
- Entrypoint routes by `OMICSAGE_MODALITY` env var or first CLI argument
- Data/reports always volume-mounted, never baked into image
- API keys injected at runtime via env vars, never stored in image
- `HDF5_USE_FILE_LOCKING=FALSE` set globally (WSL2 compatibility)
- GPU opt-in via `OMICSAGE_GPU=1` env var

---

## Current Test Count
**1472 passing, 58 skipped** (as of 2026-06-16)

---

## What Is NOT Yet Built
- Streamlit no-code UI ← NEXT
- Nextflow DSL2 orchestration
- CLI (Click-based project manager)
- Paper / benchmarking
- omicsage.io website
