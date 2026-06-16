# OmicSage — Next Session
> Written: 2026-06-16
> Phase: 9 — Streamlit UI

---

## Session Context

**Last thing completed:**
Docker containerization — Phase 8 complete.
All four modalities containerized in a single image with a modality-agnostic
entrypoint. Verified build succeeds and pipeline runs inside container with
volume-mounted data and reports.

**Files added this session:**
```
docker/Dockerfile.pipeline     ← conda-based image, all 4 modalities
docker/Dockerfile.dev          ← extends pipeline image + Jupyter + pytest
docker/entrypoint.sh           ← routes scrna/cite/multiome/spatial + jupyter/pytest/bash
docker-compose.yml             ← 4 named services + dev service
.dockerignore                  ← excludes data/, caches, .git
scripts/run_docker.sh          ← convenience wrapper with volume mounts + GPU flag
```

**Current test count:**
```bash
conda activate omicsage
python -m pytest tests/ -q --tb=short
```
Expected: 1472 passing, 58 skipped. Confirm before writing any new code.

**Docker verification commands:**
```bash
# Build (already done — use cache)
docker build -f docker/Dockerfile.pipeline -t omicsage:latest .

# Run one step to verify
./scripts/run_docker.sh scrna \
  --config config/runs/gse166635_hcc.yaml --step ingest

./scripts/run_docker.sh spatial \
  --config config/runs/kuppe_heart.yaml --step ingest
```

---

## Today's Goal

**Single deliverable: Streamlit UI — skeleton + scRNA workflow.**

A working `streamlit run ui/app.py` that lets a user:
1. Upload a dataset (10x MTX folder or H5AD file)
2. Select modality (scRNA / CITE-seq / Multiome / Spatial)
3. Configure key parameters via sliders/dropdowns (no YAML editing)
4. Run the pipeline with a live progress log
5. View the combined HTML report inline when done

Scope for this session: **scRNA modality only**. Other modalities are
placeholders. The architecture must be extensible so adding them later
is trivial.

---

## Streamlit UI Architecture

```
ui/
  app.py                  ← entry point, sidebar navigation
  pages/
    01_upload.py          ← drag-drop upload + modality selector
    02_configure.py       ← parameter sliders, generates config YAML
    03_run.py             ← pipeline runner, live log streaming
    04_report.py          ← renders combined HTML report inline
  components/
    sidebar.py            ← shared sidebar (project selector, status)
    config_builder.py     ← YAML config generator from UI inputs
    log_streamer.py       ← subprocess stdout → st.empty() live log
  state.py                ← st.session_state keys and defaults
```

**Key Streamlit decisions:**
- Use `st.session_state` for all inter-page state (config, run status)
- Pipeline runs as a `subprocess.Popen` with stdout piped to UI
- Config generated as a temp YAML file, passed to runner via `--config`
- Report rendered via `st.components.v1.html()` with the combined HTML

---

## Design Principles for the UI

- **Biologist-first**: no YAML, no terminal, no Python knowledge required
- **Sane defaults**: every parameter has a sensible default pre-filled
- **Progressive disclosure**: advanced parameters collapsed by default
- **Non-blocking**: pipeline runs in background, UI stays responsive
- **Graceful degradation**: AI features shown as optional, clearly labelled

---

## Known Issues / Watch Out For

- Streamlit reruns the entire script on every widget interaction —
  use `st.session_state` carefully to avoid resetting pipeline state
- `subprocess.Popen` stdout streaming requires `bufsize=1` and
  `universal_newlines=True` to get line-by-line output
- Large H5AD uploads: use `st.file_uploader` with `type=["h5ad","h5"]`
  and save to a temp dir, not memory
- The combined HTML report contains embedded JS — render with
  `st.components.v1.html(html, height=800, scrolling=True)`
- Do NOT import scanpy/anndata at the top of app.py —
  Streamlit reruns are slow enough without loading the full stack on
  every widget click. Import inside the run callback only.

---

## Pre-session Checklist

- [ ] Confirm baseline tests still pass (1472)
- [ ] Confirm `streamlit` is in the conda env:
      `conda activate omicsage && streamlit --version`
- [ ] If not installed: `pip install streamlit`
- [ ] Review `run_scrna_pipeline.py --help` to confirm CLI flags
      available for the UI to pass through

---

## Verify at End of Session

```bash
conda activate omicsage
streamlit run ui/app.py

# In browser: upload a small H5AD, configure, run --step ingest,
# confirm checkpoint appears and log streams live
```

---

## Memory Files to Update at End of This Session

```
.dev_memory/NEXT_SESSION.md    ← replace with next phase content
.dev_memory/CURRENT_STATUS.md  ← add Streamlit UI to completed
.dev_memory/PROGRESS.md        ← tick Streamlit phase
```
