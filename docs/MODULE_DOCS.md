# OmicSage — QC Module: Script Reference

> Location: `pipeline/modules/qc/`
> Phase: 1 — Core scRNA Pipeline
> Last updated: May 2026

This document describes every script in the QC module — what it does,
what goes in, what comes out, and how it connects to the next step.

---

## 1. `ingest.py`

**What it does**

Loads any supported single-cell dataset file and returns a clean AnnData
object with **raw integer counts guaranteed to be in `adata.X`**, regardless
of how the original file stored them.

**Why it exists**

Public datasets (e.g. GSE194122) are often distributed as processed h5ad files
where `adata.X` holds normalized values and raw counts are hidden in
`adata.layers['counts']` or `adata.raw`. Without this step, every downstream
module would need its own format-detection logic. `ingest.py` solves this once
so QC, normalization, and all other modules can safely assume raw counts in `X`.

**Supported input formats**

| Format | How to pass it |
|--------|---------------|
| `.h5ad` | Path to AnnData file |
| `.h5`   | 10x Genomics HDF5 |
| directory | 10x MEX folder (barcodes + features + matrix) |

**Key behaviours**

- Auto-detects format from file extension or directory contents
- For processed h5ad files: searches `adata.raw`, then `adata.layers` for integer counts; moves them to `adata.X` and saves the normalized matrix to `adata.layers['normalized']`
- Converts `adata.X` to sparse CSR format for memory efficiency
- Attaches provenance metadata to `adata.uns['omicsage_source']`
- Adds `adata.obs['sample']` label

**Input**

```
path  : str | Path   →  .h5ad file, .h5 file, or 10x MEX directory
```

**Output**

```
adata.X                       →  raw integer counts (CSR sparse)
adata.layers['normalized']    →  original normalized matrix (if input was processed)
adata.obs['sample']           →  sample name
adata.uns['omicsage_source']  →  format, file path, raw_layer provenance
```

**Usage**

```python
from pipeline.modules.qc.ingest import load_dataset

adata = load_dataset("data/benchmark/GSE194122_BMMC.h5ad")
adata = load_dataset("data/benchmark/sample.h5")
adata = load_dataset("data/benchmark/HCC1/")   # MTX directory
```

**Connects to**: `qc.py` — pass the returned AnnData directly to `run_qc()`

---

## 2. `qc.py`

**What it does**

Performs quality control on a raw-count AnnData object. Computes per-cell
metrics, detects doublets, applies configurable filters, and returns a clean
filtered AnnData plus a summary metrics dictionary.

**Why it exists**

Low-quality cells (dead cells, empty droplets, multiplets) distort clustering,
differential expression, and every downstream analysis. This step removes them
using standard and well-validated criteria before any biological analysis begins.

**Steps performed (in order)**

1. `var_names_make_unique()` — handles CITE-seq files with mixed RNA + ADT features
2. Mitochondrial gene detection — auto-detects `MT-` (human) or `mt-` (mouse) prefix; falls back to `adata.var['gene_ids']` if var_names don't contain symbols
3. Per-cell QC metrics via `sc.pp.calculate_qc_metrics()`:
   - `n_genes_by_counts` — genes detected per cell
   - `total_counts` — total UMI count per cell
   - `pct_counts_mt` — mitochondrial read percentage
4. Doublet detection via Scrublet — adds `obs['doublet_score']` and `obs['predicted_doublet']`; fails gracefully if Scrublet errors
5. Filter application — removes cells failing any threshold
6. Optional HTML report generation via `qc_report.py`

**Default thresholds**

| Parameter | Default | What it removes |
|-----------|---------|----------------|
| `min_genes` | 200 | Empty droplets / dead cells |
| `max_genes` | 2500 | Likely multiplets |
| `max_mt_pct` | 5% | Lysed / dying cells |
| `remove_doublets` | True | Scrublet-detected doublets |

All thresholds are configurable per-dataset.

**Input**

```
adata  : AnnData   →  raw counts in adata.X (output of ingest.py)
```

**Output**

```
adata_filtered  : AnnData   →  cells passing all QC filters; QC metrics in .obs
metrics         : dict      →  cell counts, per-filter removal counts, medians, thresholds
```

**Usage**

```python
from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc

adata = load_dataset("data/benchmark/GSE194122_BMMC.h5ad")

# Basic usage
adata_filtered, metrics = run_qc(adata)

# Custom thresholds + HTML report
adata_filtered, metrics = run_qc(
    adata,
    min_genes=300,
    max_genes=5000,
    max_mt_pct=15.0,
    generate_report=True,
    report_path="reports/qc_GSE194122.html",
    sample_name="GSE194122_BMMC",
)
```

**Connects to**: `normalization.py` (Phase 1, next step) — pass `adata_filtered`

---

## 3. `qc_report.py`

**What it does**

Generates a self-contained HTML QC report after filtering. The report
is a single `.html` file with all plots embedded as base64 PNG images —
no external dependencies, no internet required, opens in any browser.

**Why it exists**

OmicSage's core promise is that every analysis step produces a readable
output — not just a filtered AnnData that only a bioinformatician can
inspect. The QC report lets a biologist (or a PI reviewing the analysis)
understand exactly how many cells were removed and why, without looking
at code.

**Report contents**

| Section | What it shows |
|---------|--------------|
| Summary cards | Cells in, cells kept, cells removed, pass rate, MT genes found, doublets removed |
| Violin plots (before vs after) | Genes per cell, total UMI, MT% — with threshold lines |
| UMI vs genes scatter | Each dot = one cell, coloured by MT% |
| Doublet score histogram | Scrublet score distribution with threshold marker |
| Distribution medians table | Median genes, UMI, MT% pre-QC |
| Filter thresholds table | Every parameter used + cells removed by that filter |
| Ground-truth validation | MT% correlation vs `GEX_pct_counts_mt` (only if column present) |

**Input**

```
adata_raw      : AnnData  →  pre-filter AnnData (after metrics computed)
adata_filtered : AnnData  →  post-filter AnnData
metrics        : dict     →  output of run_qc()
output_path    : str      →  where to write the .html file
sample_name    : str      →  label shown in the report header
```

**Output**

```
reports/qc_report.html   →  self-contained HTML file (~1-3 MB)
```

**Usage**

Called automatically by `run_qc()` when `generate_report=True`.
Can also be called directly:

```python
from pipeline.modules.qc.qc_report import generate_qc_report

generate_qc_report(
    adata_raw=adata_before_filter,
    adata_filtered=adata_after_filter,
    metrics=metrics,
    output_path="reports/qc_GSE194122.html",
    sample_name="GSE194122_BMMC",
)
```

**Dependencies**: `matplotlib` only — no Plotly, no Bokeh, no JS frameworks.

---

## 4. `data_report.py`

**What it does**

Generates a data intake HTML report summarising the contents of a raw
h5ad file *before any QC is applied*. Designed to be run immediately
after downloading a new dataset to understand its structure.

**Why it exists**

Before running any pipeline step, you need to know: how many cells, how
many genes, what metadata columns exist, what layers are present, and
whether the file needs special handling. This report answers all of
those questions in one browser tab.

**Important**: Uses `backed='r'` mode to read h5ad files — only metadata
is loaded into memory. The full count matrix is never read. This means
the report works on files of any size (including 2.9 GB multiome files)
without memory issues.

**Input**

```
--input   : path to .h5ad file
--geo     : GEO accession (e.g. GSE194122) shown in report header
--output  : path for HTML output (e.g. reports/data_intake.html)
```

**Output**

```
reports/data_intake.html   →  self-contained HTML file
```

**Usage**

```bash
python pipeline/modules/qc/data_report.py \
  --input  data/benchmark/GSE194122_BMMC.h5ad \
  --geo    GSE194122 \
  --output reports/data_intake_BMMC.html
```

**Connects to**: Read the report → decide QC thresholds → run `qc.py`

---

## Module Data Flow

```
Raw file (.h5ad / .h5 / MTX dir)
        │
        ▼
   ingest.py          → adata.X = raw counts
        │
        ▼
     qc.py            → adata_filtered + metrics dict
        │                       │
        │                       ▼
        │               qc_report.py  → reports/qc_report.html
        │
        ▼
  normalization.py    ← NEXT STEP (Phase 1)
```

---

## Running the Full QC Step

```bash
cd ~/OmicSage
conda activate omicsage

python3 -c "
from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc

adata = load_dataset('data/benchmark/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad')
adata_filtered, metrics = run_qc(
    adata,
    generate_report=True,
    report_path='reports/qc_BMMC.html',
    sample_name='GSE194122_BMMC',
)
print('Cells before:', metrics['n_cells_input'])
print('Cells after: ', metrics['n_cells_output'])
print('Removed:     ', metrics['n_cells_removed'])
"
```

---

## Tests

| Test file | What it covers |
|-----------|---------------|
| `tests/test_phase0_structure.py` | Repo structure, imports, config schema |
| `tests/test_ingest.py` | Format detection, raw count extraction, all three loaders |
| `tests/test_qc.py` | MT detection, metric computation, filtering, Scrublet, ground-truth validation |

Run all QC-related tests:

```bash
python -m pytest tests/test_phase0_structure.py tests/test_ingest.py tests/test_qc.py -v
```
