# OmicSage — Next Session
> Written: 2026-06-11
> Phase: 8 — Docker Containerization

---

## Session Context

**Last thing completed:**
- Phase 7 fully complete (Sessions 1–5 + Session A + Session B extensions)
- Session A: `spatial_impute.py` + `spatial_impute_report.py` + `test_spatial_impute.py`
  + runner + combined report + kuppe config updated. Multiple bug fixes:
  OOM fix (Tangram clusters mode), `project_genes` API fix, numpy array
  storage, log-normalisation in validation scatter, library_key multi-sample fix.
- Session B: `spatial_ingest.py` Visium HD + Xenium loaders (spatialdata-io),
  `test_spatial_ingest_formats.py` (27 tests), `visium_hd_mouse_brain.yaml`,
  `xenium_breast.yaml`, `scripts/download_spatial_benchmark.py`.
- Docs: `SPATIAL_MASTER_PROMPT.md` v1.1, `SPATIAL_MODULE_DOCS.md` fully updated.
- Tests fixed: `test_spatial_ingest_qc.py` (5 stub tests updated),
  `test_spatial_impute.py` (mock project_genes return value fixed).

**Current test count:**
```bash
conda activate omicsage
python -m pytest tests/ -q --tb=short
```
Confirm baseline before writing any new code. Expected: ~1476 passing.

**Files modified last session block:**
```
pipeline/modules/spatial/spatial_impute.py         ← NEW
reports/templates/spatial/spatial_impute_report.py ← NEW
reports/spatial_combined_report.py                 ← Impute tab added
run_spatial_pipeline.py                            ← impute step added
tests/test_spatial_impute.py                       ← NEW
tests/test_spatial_ingest_qc.py                    ← stub tests updated
tests/test_spatial_ingest_formats.py               ← NEW
pipeline/modules/spatial/spatial_ingest.py         ← Visium HD + Xenium loaders
config/runs/kuppe_heart.yaml                       ← impute block added
config/runs/visium_hd_mouse_brain.yaml             ← NEW
config/runs/xenium_breast.yaml                     ← NEW
scripts/download_spatial_benchmark.py              ← NEW
ai/manual_review/SPATIAL_MASTER_PROMPT.md          ← v1.1
SPATIAL_MODULE_DOCS.md                             ← fully updated
```

---

## Today's Goal

**Single deliverable: Docker containerization of the OmicSage pipeline.**

A working `docker-compose up` that runs the full spatial pipeline end-to-end
on the Kuppe benchmark dataset inside a container, with results written to the
host filesystem.

---

## Deliverables

```
docker/
  Dockerfile.pipeline       ← Python analysis environment (conda-based)
  Dockerfile.dev            ← development image with Jupyter + hot-reload
docker-compose.yml          ← orchestrates pipeline + volume mounts
.dockerignore               ← exclude data/, checkpoints/, __pycache__ etc.
scripts/run_docker.sh       ← convenience wrapper: docker run with correct mounts
```

Optional but useful this session:
```
docker/entrypoint.sh        ← flexible entrypoint: run pipeline step or shell
```

---

## Design Decisions to Make Before Coding

**1. Base image**
Use `continuumio/miniconda3` as the base — matches your WSL2 conda environment
exactly. Alternative is `python:3.11-slim` + pip-only, but conda is safer given
the heavy bioinformatics dependency stack (scanpy, squidpy, tangram-sc, etc.).

**2. Dependency installation strategy**
Two-stage build:
- Stage 1: install `environment.yml` into conda env `omicsage`
  (`conda env create -f environment.yml`)
- Stage 2: copy source code
This keeps the layer cache efficient — a code change doesn't invalidate
the 3 GB dependency layer.

**3. Volume mounts (critical)**
```
data/benchmark/   → read-only input
data/processed/   → read-write output (checkpoints)
reports/          → read-write output (HTML reports)
config/           → read-only config
pipeline/         → read-write during dev (hot-reload)
```
Never bake data into the image.

**4. Entrypoint**
```bash
conda run -n omicsage python run_spatial_pipeline.py \
  --config config/runs/kuppe_heart.yaml \
  "$@"          # pass --step, --from-step, --force etc. through
```

**5. GPU support**
Add `--gpus all` flag to `docker run` in the convenience script,
guarded by `OMICSAGE_GPU=1` env var. Default is CPU-only.

---

## Dockerfile.pipeline skeleton

```dockerfile
FROM continuumio/miniconda3:latest

# System deps needed by some bio packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ libhdf5-dev git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Layer 1: dependencies (cached unless environment.yml changes)
COPY environment.yml .
RUN conda env create -f environment.yml \
    && conda clean -afy

# Layer 2: source code
COPY . .

# Ensure conda env is on PATH
ENV PATH /opt/conda/envs/omicsage/bin:$PATH

ENTRYPOINT ["python", "run_spatial_pipeline.py"]
CMD ["--help"]
```

---

## docker-compose.yml skeleton

```yaml
version: "3.9"
services:
  pipeline:
    build:
      context: .
      dockerfile: docker/Dockerfile.pipeline
    volumes:
      - ./data:/app/data
      - ./reports:/app/reports
      - ./config:/app/config
      - ./pipeline:/app/pipeline     # hot-reload during dev
    environment:
      - OMICSAGE_GPU=${OMICSAGE_GPU:-0}
    command: ["--config", "config/runs/kuppe_heart.yaml"]
```

---

## Pre-session checklist

Before starting:
- [ ] Run baseline test count and record it here
- [ ] Confirm Docker Desktop is running in WSL2
      (`docker info` should return without error)
- [ ] Confirm `environment.yml` is up to date — check that these are present:
      `tangram-sc`, `spatialdata-io`, `squidpy`, `scvi-tools`, `cell2location`
      (they should be, but verify before baking into the image)
- [ ] Note the size of the conda env:
      `du -sh ~/miniconda3/envs/omicsage/`
      This is roughly the Docker image size — plan for a 4–8 GB image

---

## Known Issues / Watch Out For

- `continuumio/miniconda3` base images can be slow to pull (~150 MB).
  Build once, iterate on the source layer.
- `conda env create` inside Docker can hit memory limits on WSL2.
  If it OOMs, add `--memory 8g` to the Docker build command.
- Some packages (cell2location, tangram-sc) compile C extensions — the
  `build-essential` apt layer is required.
- `hdf5` headers (`libhdf5-dev`) are needed for `h5py` native compilation.
- On WSL2, Docker Desktop must have "Use WSL 2 based engine" enabled
  (Docker Desktop → Settings → General).
- `.dockerignore` must exclude `data/` to keep the build context small —
  benchmark datasets are several GB and must never be baked into the image.
- `conda run -n omicsage` vs activating the env: inside Docker, activation
  doesn't work in non-interactive shells. Use the full path
  `/opt/conda/envs/omicsage/bin/python` or `ENV PATH` in the Dockerfile.

---

## Verify at End of Session

```bash
# Build the image
docker build -f docker/Dockerfile.pipeline -t omicsage:latest .

# Run one pipeline step on Kuppe (fast — ingest only)
docker run --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/reports:/app/reports \
  -v $(pwd)/config:/app/config \
  omicsage:latest \
  --config config/runs/kuppe_heart.yaml --step ingest

# Confirm checkpoint written to host filesystem
ls data/processed/kuppe_heart/01_ingested.h5ad
```

---

## Memory Files to Update at End of This Session

```
.dev_memory/NEXT_SESSION.md    ← replace with Streamlit phase content
.dev_memory/CURRENT_STATUS.md  ← add Docker to completed infrastructure
.dev_memory/PROGRESS.md        ← tick Docker phase
```
