# 🧬 OmicSage

> **AI-Assisted Single-Cell Multi-Omics Analysis Platform**
>
> End-to-end pipeline · Automated reports · AI interpretation · No API key required

[![CI](https://github.com/fshokor/OmicSage/actions/workflows/ci.yml/badge.svg)](https://github.com/fshokor/OmicSage/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Status: Alpha](https://img.shields.io/badge/status-alpha--dev-orange.svg)]()

---

## What Is OmicSage?

OmicSage is an open-source platform that covers the full single-cell multi-omics stack in one place:

| Modality | Status |
|----------|--------|
| scRNA-seq | 🔧 Phase 1 — In progress |
| scATAC-seq | 📅 Phase 4 — Planned |
| Spatial transcriptomics | 📅 Phase 5 — Planned |
| Multiome (RNA + ATAC) | 📅 Phase 6 — Planned |

**Two users. One tool.**
- 🔬 **Biologists**: no-code Streamlit interface, guided workflow, biological interpretation in plain language
- 💻 **Bioinformaticians**: project templates, multi-project management, automated reports, no repeated code

---

## Architecture

```
┌─────────────────────────────────────────────────────┐
│  LAYER 3: AI Intelligence (OPTIONAL)                │
│  BioChatter + LLM → threshold suggestions,          │
│  cluster interpretation, PubMed RAG, narratives     │
├─────────────────────────────────────────────────────┤
│  LAYER 2: Report Engine (ALWAYS ON)                 │
│  Quarto + python-pptx → HTML/PDF reports +          │
│  PowerPoint slides after every analysis step        │
├─────────────────────────────────────────────────────┤
│  LAYER 1: Core Pipeline (ALWAYS ON)                 │
│  Nextflow DSL2 → QC → normalize → integrate →       │
│  cluster → annotate → downstream analysis           │
└─────────────────────────────────────────────────────┘
```

**Key principle**: `ai_features: false` in your config runs the full pipeline without any API key or internet connection.

---

## Quick Start

### Requirements
- Docker Desktop + WSL2 (Windows) **or** Docker (Linux/macOS)
- Nextflow ≥ 23.04

### 1. Clone the repo
```bash
git clone https://github.com/fshokor/OmicSage.git
cd OmicSage
```

### 2. Build the Docker images
```bash
docker build -f docker/Dockerfile.python -t omicsage/python:latest .
# Optional: R-based steps (Seurat, SingleR, DESeq2)
docker build -f docker/Dockerfile.r -t omicsage/r:latest .
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
python cli/omicsage.py run my_analysis/ --profile docker
```

### Web UI (biologists)
```bash
docker compose up streamlit
# Open http://localhost:8501
```

---

## AI Layer

OmicSage uses [BioChatter](https://github.com/biocypher/biochatter) for AI features. Set your provider in `config.yaml`:

```yaml
ai:
  enabled: true
  provider: ollama         # fully local, no API key
  ollama_model: llama3.2
```

Or with Claude:
```bash
export ANTHROPIC_API_KEY=sk-ant-...
# config.yaml: provider: claude, model: claude-opus-4-5
```

**All AI calls are audit-logged** to `logs/llm/` in JSONL format.

---

## Roadmap

| Phase | Scope | Target Week |
|-------|-------|------------|
| 0 ✅ | Foundation — repo, Docker, CI/CD | 1-2 |
| 1 🔧 | Core scRNA-seq pipeline | 2-6 |
| 2 | Report engine (Quarto + PowerPoint) | 6-9 |
| 3 | AI layer (BioChatter integration) | 9-13 |
| 4 | scATAC-seq module | 13-16 |
| 5 | Spatial transcriptomics | 16-19 |
| 6 | Multiome integration | 19-22 |
| 7 | Streamlit web UI | 22-25 |
| 8 | Benchmark + paper | 25-30 |

---

## Citation

OmicSage is under active development. A preprint will be posted on bioRxiv in 2026. If you use OmicSage in your research, please check back for citation details.

---

## License

MIT — see [LICENSE](LICENSE)
