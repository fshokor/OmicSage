# 🧬 OmicSage

> **Single-Cell Multi-Omics Analysis Platform**
>
> End-to-end pipeline · Automated reports · Guided interpretation · Full Python

[![CI](https://github.com/fshokor/OmicSage/actions/workflows/ci.yml/badge.svg)](https://github.com/fshokor/OmicSage/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/scanpy-1.10-blue.svg)](https://scanpy.readthedocs.io/)
[![Status: Alpha](https://img.shields.io/badge/status-alpha--dev-orange.svg)]()

---

## What Is OmicSage?

OmicSage is an open-source platform that covers the full single-cell multi-omics stack in one place:

| Modality | Status |
|----------|--------|
| scRNA-seq | ✅ Phase 1 —  Finished |
| Cite-seq | ✅ Phase 4 — Finished |
| Multiome (RNA + ATAC) | 🔧 Phase 5 — In progress |
| Spatial transcriptomics | 📅 Phase 6 — Planned |

**Two users. One tool.**
- 🔬 **Biologists**: no-code Streamlit interface, guided workflow, biological interpretation in plain language
- 💻 **Bioinformaticians**: project templates, multi-project management, automated reports, no repeated code

---

## Architecture

```
┌─────────────────────────────────────────────────────┐
│  LAYER 3: Review & Interpretation (MANUAL)          │
│  Structured outputs + QC flags guide the analyst    │
│  through threshold decisions and cluster review     │
├─────────────────────────────────────────────────────┤
│  LAYER 2: Report Engine (ALWAYS ON)                 │
│  Per-step HTML reports + combined tabbed report     │
│  → 00_combined_report.html after every pipeline run │
├─────────────────────────────────────────────────────┤
│  LAYER 1: Core Pipeline (ALWAYS ON)                 │
│  Python (Scanpy/scverse) → QC → normalize →         │
│  integrate → cluster → annotate → downstream        │
└─────────────────────────────────────────────────────┘
```

**Key principle**: every step produces structured outputs and a human-readable report so the analyst stays in control of all biological decisions.

---

## Quick Start

### Requirements
- Python 3.11
- Conda (Miniconda or Anaconda)

### 1. Clone the repo
```bash
git clone https://github.com/fshokor/OmicSage.git
cd OmicSage
```

### 2. Create and activate the conda environment
```bash
conda env create -f environment.yml
conda activate omicsage
```

### 3. Create your first project
```bash
python cli/omicsage.py create-project my_analysis --modality scrna
```

### 4. Edit the config
```bash
nano my_analysis/config.yaml   # set input.scrna_path to your data
```

### 5. Run the pipeline
```bash
python cli/omicsage.py run my_analysis/
```

### Web UI (biologists)
```bash
conda activate omicsage
streamlit run ui/app.py
# Open http://localhost:8501
```

---

## Running Tests

```bash
conda activate omicsage
python -m pytest tests/ -v
```

> Always use `python -m pytest`, not bare `pytest`, to ensure the correct Python environment is used.

---

## Interpretation Layer

OmicSage is designed around **manual review by the analyst**. Rather than delegating biological decisions to an LLM, every pipeline step produces structured outputs and a human-readable HTML report that guides you through the results:

- **QC step** — per-sample MAD-based thresholds surfaced with flags; analyst confirms or adjusts
- **Clustering step** — UMAP + marker gene tables per cluster for manual label assignment
- **Annotation step** — cell type predictions from [CelltypistML](https://github.com/Teichlab/celltypist) (Python) alongside marker evidence; analyst reviews and overrides
- **DEG step** — ranked gene tables with volcano plots; analyst interprets biological significance
- **Combined report** — `00_combined_report.html` collects all steps in one tabbed view after every run

This keeps the analyst in control of all biological decisions while eliminating boilerplate code and repetitive plotting.

---

## Roadmap

| Phase | Scope | Target Week |
|-------|-------|------------|
| 0 ✅ | Foundation — repo, CI/CD | 1-2 |
| 1 ✅ | Core scRNA-seq pipeline (annotation in progress) | 2-6 |
| 2 ✅ | Report engine — combined tabbed HTML report | 6-9 |
| 3 ✅ | Manual review layer — structured QC flags, cluster review, CelltypistML annotation | 9-13 |
| 4 ✅ | Cite-seq module | 13-16 |
| 5 🔧 | Multiome integration | 16-19 |
| 6 | Spatial transcriptomics | 19-22 |
| 7 | Streamlit web UI | 22-25 |
| 8 | Benchmark + paper | 25-30 |

---

## Citation

OmicSage is under active development. A preprint will be posted on bioRxiv in 2026. If you use OmicSage in your research, please check back for citation details.

---

## License

MIT — see [LICENSE](LICENSE)
