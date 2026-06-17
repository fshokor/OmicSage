# 🧬 OmicSage

> **Open-Source Single-Cell Multi-Omics Analysis Platform**
>
> End-to-end pipeline · Automated HTML reports · No-code Streamlit UI · HPC-ready via Nextflow

[![CI](https://github.com/fshokor/OmicSage/actions/workflows/ci.yml/badge.svg)](https://github.com/fshokor/OmicSage/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen.svg)](https://www.nextflow.io/)
[![Status: Alpha](https://img.shields.io/badge/status-alpha--dev-orange.svg)]()

---

<!-- SCREENSHOT: Streamlit UI — Dataset page -->
<!-- Replace with: docs/images/ui_dataset.png -->
![OmicSage Streamlit UI](docs/images/placeholder_ui.png)

---

## What Is OmicSage?

OmicSage is an open-source platform that covers the full single-cell multi-omics analysis stack in one place — from raw count matrices to publication-ready reports.

**No existing tool delivers all of this together:**

| Feature | OmicSage | nf-core/scdownstream | CellWhisperer | SCAPE |
|---------|----------|---------------------|---------------|-------|
| scRNA-seq pipeline | ✅ | ✅ | ❌ | ✅ |
| CITE-seq pipeline | ✅ | ❌ | ❌ | partial |
| Multiome (RNA+ATAC) | ✅ | ❌ | ❌ | partial |
| Spatial transcriptomics | ✅ | ❌ | ❌ | ❌ |
| Automated HTML reports | ✅ | ❌ | ❌ | ❌ |
| No-code UI | ✅ | ❌ | ✅ | ❌ |
| HPC/cloud (Nextflow) | ✅ | ✅ | ❌ | ❌ |
| Works without API key | ✅ | ✅ | ❌ | ❌ |

**Two users. One tool.**
- 🔬 **Biologists** — no-code Streamlit interface, step-by-step workflow, biological interpretation in plain language
- 💻 **Bioinformaticians** — config-driven runners, multi-project management, automated reports, HPC submission via Nextflow

---

## Supported Modalities

| Modality | Steps | Runner |
|----------|-------|--------|
| **scRNA-seq** | QC → normalize → reduce → cluster → annotate → DEG → GSEA → Harmony → pseudobulk | `run_scrna_pipeline.py` |
| **CITE-seq** | ADT normalize → doublets → reduce → harmony → annotate → integration → DEG → GSEA → correlation → epitope | `run_cite_pipeline.py` |
| **Multiome (RNA+ATAC)** | ATAC QC → LSI reduce → annotate → MultiVI integration → DEG → GRN | `run_multiome_pipeline.py` |
| **Spatial** | Ingest → QC → reduce → cluster → deconvolve → downstream → impute | `run_spatial_pipeline.py` |

---

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  INTERFACE LAYER                                            │
│  Streamlit UI (biologists) · Python CLI · Nextflow (HPC)   │
├─────────────────────────────────────────────────────────────┤
│  REPORT ENGINE (always on)                                  │
│  Per-step HTML reports → 00_combined_report.html            │
│  Generated automatically after every pipeline run           │
├─────────────────────────────────────────────────────────────┤
│  CORE PIPELINE (always on)                                  │
│  Scanpy · Muon · squidpy · scvi-tools · decoupler           │
│  QC → normalize → integrate → cluster → annotate →         │
│  DEG → GSEA → GRN → spatial downstream                     │
├─────────────────────────────────────────────────────────────┤
│  INFRASTRUCTURE                                             │
│  Docker (omicsage:latest) · Nextflow DSL2 · GitHub Actions  │
└─────────────────────────────────────────────────────────────┘
```

---

## Quick Start

### Requirements
- Python 3.11
- Conda (Miniconda or Anaconda)
- Docker Desktop (optional — for containerized runs)
- Java 21 + Nextflow (optional — for HPC/cloud)

### Installation

```bash
git clone https://github.com/fshokor/OmicSage.git
cd OmicSage
conda env create -f environment.yml
conda activate omicsage
```

---

### Path A — Streamlit UI (biologists, no coding required)

```bash
conda activate omicsage
streamlit run ui/app.py
# Open http://localhost:8501
```

<!-- SCREENSHOT: Streamlit Configure page -->
<!-- Replace with: docs/images/ui_configure.png -->
![Configure Page](docs/images/placeholder_configure.png)

1. **Dataset** — enter your dataset name, select modality, set the path to your data
2. **Configure** — toggle steps on/off, adjust parameters via sliders and dropdowns
3. **Run** — click ▶ Run Pipeline, watch live log output
4. **Report** — view and download the combined HTML report

<!-- SCREENSHOT: Streamlit Run page with live log -->
<!-- Replace with: docs/images/ui_run.png -->
![Run Page](docs/images/placeholder_run.png)

---

### Path B — Python CLI (bioinformaticians)

```bash
conda activate omicsage

# Run the full scRNA pipeline
python run_scrna_pipeline.py --config config/runs/GSE166635.yaml

# Run a single step
python run_scrna_pipeline.py --config config/runs/GSE166635.yaml --step qc

# Run a range of steps
python run_scrna_pipeline.py --config config/runs/GSE166635.yaml \
  --from-step normalize --to-step annotate

# Re-run from a step (force overwrite checkpoints)
python run_scrna_pipeline.py --config config/runs/GSE166635.yaml \
  --from-step cluster --force

# Same pattern for other modalities
python run_cite_pipeline.py     --config config/runs/GSE194122_cite.yaml
python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml
python run_spatial_pipeline.py  --config config/runs/kuppe_heart.yaml
```

Each step writes a checkpoint `.h5ad` file to `processed_dir` and an HTML report to `reports_dir`.

---

### Path C — Nextflow (HPC / cloud)

```bash
# Install Nextflow (requires Java 21)
curl -s https://get.nextflow.io | bash

# Run full pipeline locally with Docker
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile local

# Resume after a crash (skips completed steps)
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile local \
  -resume

# Run on SLURM (HPC)
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile slurm

# Run with GPU (totalVI, Tangram, MultiVI)
nextflow run main.nf \
  --config config/runs/GSE194122_cite.yaml \
  --modality cite \
  -profile local \
  --gpu true
```

---

### Path D — Docker

```bash
# Build the image
docker build -f docker/Dockerfile.pipeline -t omicsage:latest .

# Run any modality
docker run --rm \
  -v $(pwd):/app \
  --entrypoint '' \
  omicsage:latest \
  /opt/conda/envs/omicsage/bin/python /app/run_scrna_pipeline.py \
  --config /app/config/runs/GSE166635.yaml

# With GPU
docker run --rm --gpus all \
  -e OMICSAGE_GPU=1 \
  -v $(pwd):/app \
  --entrypoint '' \
  omicsage:latest \
  /opt/conda/envs/omicsage/bin/python /app/run_cite_pipeline.py \
  --config /app/config/runs/GSE194122_cite.yaml
```

---

## Configuration

Every analysis is driven by a YAML config file. Example configs are in `config/runs/`:

```
config/runs/
  GSE166635.yaml              ← HCC scRNA-seq (Wang et al. 2025)
  GSE194122_cite.yaml         ← BMMC CITE-seq (NeurIPS 2021)
  GSE194122_atac_multiome.yaml← BMMC Multiome (NeurIPS 2021)
  kuppe_heart.yaml            ← Human heart Visium (Kuppe et al. 2022)
```

Minimal scRNA config:

```yaml
dataset:
  id: my_dataset
  name: "My scRNA-seq dataset"

paths:
  raw_input: data/raw/my_dataset
  processed_dir: data/processed/my_dataset
  reports_dir: reports/my_dataset

steps:
  qc:
    enabled: true
    params:
      min_genes: 200
      max_genes: 6000
      max_mt_pct: 20.0
      remove_doublets: true
  normalize:
    enabled: true
  # ... add steps as needed
```

---

## Outputs

After a pipeline run, you get:

```
data/processed/my_dataset/
  01_qc.h5ad
  02_normalized.h5ad
  03_reduced.h5ad
  04_clustered.h5ad
  05_annotated.h5ad
  06_deg.h5ad
  07_gsea.h5ad
  ...

reports/my_dataset/
  01_qc_report.html
  02_normalize_report.html
  03_reduce_report.html
  04_cluster_report.html
  05_annotate_report.html
  06_deg_report.html
  07_gsea_report.html
  00_combined_report.html     ← all tabs in one file
```

<!-- SCREENSHOT: Combined HTML report — QC tab -->
<!-- Replace with: docs/images/report_qc.png -->
![Combined Report QC Tab](docs/images/placeholder_report_qc.png)

<!-- SCREENSHOT: Combined HTML report — UMAP/cluster tab -->
<!-- Replace with: docs/images/report_cluster.png -->
![Combined Report Cluster Tab](docs/images/placeholder_report_cluster.png)

---

## Running Tests

```bash
conda activate omicsage
python -m pytest tests/ -q --tb=short
# Expected: 1472 passing, 58 skipped
```

---

## What's Built

| Component | Status |
|-----------|--------|
| scRNA-seq pipeline (10 steps) | ✅ Complete |
| CITE-seq pipeline (10 steps) | ✅ Complete |
| Multiome pipeline (6 steps) | ✅ Complete |
| Spatial pipeline (7 steps) | ✅ Complete |
| Per-step HTML reports | ✅ Complete |
| Combined tabbed report | ✅ Complete |
| AI interpretation layer | ✅ Built (optional, off by default) |
| Streamlit UI | ✅ Complete |
| Docker containerization | ✅ Complete |
| Nextflow DSL2 (all 4 modalities) | ✅ Complete |
| Tutorial notebook | 🔧 In progress |
| Google Colab notebook | 🔧 In progress |
| CLI (Click-based) | ⏳ Planned |
| omicsage.io website | ⏳ Planned |
| Paper / preprint | ⏳ Planned |

---

## Benchmark Datasets

OmicSage is validated on three published datasets:

| Dataset | Modality | Reference |
|---------|----------|-----------|
| GSE166635 | scRNA-seq | Wang et al. 2025, npj Precision Oncology (HCC) |
| GSE194122 | CITE-seq + Multiome | NeurIPS 2021 Open Problems benchmark (BMMC) |
| Kuppe heart | Spatial (Visium) | Kuppe et al. 2022, Nature (myocardial infarction) |

---

## Citation

OmicSage is under active development. A preprint will be posted on bioRxiv in 2026.
If you use OmicSage in your research, please check back for citation details.

---

## License

MIT — see [LICENSE](LICENSE)
