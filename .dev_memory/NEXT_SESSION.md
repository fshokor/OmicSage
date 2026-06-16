# OmicSage — Next Session
> Written: 2026-06-16
> Phase: 10 — Nextflow DSL2 Orchestration

---

## Session Context

**Last thing completed:**
Streamlit UI — Phase 9 complete.
Full 4-page UI built and debugged:
- Page 1: Dataset (modality, name, ID, organism, data path, RNA ref path)
- Page 2: Configure (step selector, per-step parameter widgets for all 4 modalities)
- Page 3: Run (step range, --force, live log streaming via background thread)
- Page 4: Report (embedded HTML via components.v1.html)
- Save config YAML to config/runs/<id>.yaml
- Load existing config from disk
- Project history (.omicsage_history.json) with sidebar quick-reload

**Docker status:**
Phase 8 complete. All 4 modalities containerized in a single image.
- docker/Dockerfile.pipeline  — conda-based, all modalities
- docker/Dockerfile.dev       — extends pipeline + Jupyter + pytest
- docker/entrypoint.sh        — routes modality + jupyter/pytest/bash
- docker-compose.yml          — 4 named services + dev service
- scripts/run_docker.sh       — convenience wrapper

**Current test count:**
```bash
conda activate omicsage
python -m pytest tests/ -q --tb=short
```
Expected: 1472 passing, 58 skipped. Confirm before writing any new code.

**UI files (current state):**
```
ui/
  app.py
  state.py
  defaults.py
  config_builder.py
  config_io.py          ← new this session (save/load YAML)
  history.py            ← new this session (project history)
  __init__.py
  components/
    __init__.py
    sidebar.py          ← updated this session (history panel)
  _pages/
    __init__.py
    p1_dataset.py       ← updated this session (load config + history)
    p2_configure.py     ← updated this session (save config button)
    p3_run.py           ← updated this session (history recording)
    p4_report.py
```

---

## Today's Goal

**Single deliverable: Nextflow DSL2 orchestration layer.**

Wrap the existing Python pipeline runners in Nextflow DSL2 so that:
1. Each pipeline step = one Nextflow process
2. Steps are chained via channel outputs (checkpoint h5ad files)
3. The full pipeline can be launched with:
   `nextflow run main.nf --config config/runs/GSE166635.yaml`
4. Steps skip automatically if checkpoint exists (Nextflow -resume)
5. Docker container used for execution (already built)

---

## Nextflow Architecture Plan

```
main.nf                     ← entry point, routes to correct workflow
nextflow.config             ← executor, container, resource settings
workflows/
  scrna.nf                  ← scRNA-seq workflow (chain of processes)
  cite.nf                   ← CITE-seq workflow
  multiome.nf               ← Multiome workflow
  spatial.nf                ← Spatial workflow
modules/
  scrna/
    qc.nf                   ← process QC { ... }
    normalize.nf
    reduce.nf
    cluster.nf
    annotate.nf
    deg.nf
    gsea.nf
    harmony.nf
    pseudobulk.nf
  cite/
    normalize_adt.nf
    ... (one file per step)
  multiome/
    atac_qc.nf
    ...
  spatial/
    ingest.nf
    ...
```

**Key design decisions to confirm at session start:**
1. Each Nextflow process calls the existing Python runner with `--step <name>`
   (reuse all existing pipeline code, no rewriting)
2. Checkpoint h5ad files are the channel values passed between processes
3. `--resume` uses Nextflow's built-in work dir caching
4. Container: `omicsage:latest` (already built in Phase 8)
5. Start with scRNA-seq only — other modalities follow same pattern

---

## Pre-session Checklist

- [ ] Confirm baseline tests still pass (1472)
- [ ] Confirm Nextflow is installed:
      `nextflow -version`
      If not: `curl -s https://get.nextflow.io | bash`
- [ ] Confirm Docker image exists:
      `docker images | grep omicsage`
- [ ] Review run_scrna_pipeline.py --help to confirm all --step names

---

## Verify at End of Session

```bash
# Run the scRNA pipeline through Nextflow
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  --step qc

# Check Nextflow work dir
ls work/

# Resume from normalize
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -resume
```

---

## Memory Files to Update at End of This Session

```
.dev_memory/NEXT_SESSION.md    ← replace with Phase 11 content
.dev_memory/CURRENT_STATUS.md  ← add Nextflow layer to completed
.dev_memory/PROGRESS.md        ← tick Nextflow phase
```
