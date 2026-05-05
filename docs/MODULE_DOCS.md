# OmicSage — QC Module: Script Reference

> Location: `pipeline/modules/qc/`
> Phase: 1 — Core scRNA Pipeline
> Last updated: 2026-05-05

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

Performs quality control on a raw-count AnnData object. Automatically detects
whether the input is plain RNA, CITE-seq (RNA + ADT), or Multiome (RNA + ATAC).
Computes all QC metrics on RNA features only, applies cell filters, and returns
a MuData object containing one AnnData per modality with low-quality cells removed.

**Why it exists**

Low-quality cells (dead cells, empty droplets, multiplets) distort clustering,
differential expression, and every downstream analysis. This step removes them
using standard and well-validated criteria before any biological analysis begins.

For multi-modal data (CITE-seq, Multiome), running QC on the full mixed feature
matrix produces biologically meaningless metrics — ADT protein counts inflate
`total_counts`, and ATAC peaks inflate `n_genes_by_counts`. `qc.py` isolates
the RNA layer for all metric computation and filtering, then propagates the
resulting cell keep-mask to all other modalities so no features are lost.

**Steps performed (in order)**

1. Modality detection — inspects `adata.var['feature_types']`; returns `"rna"`, `"cite"`, or `"multiome"`
2. Modality splitting — separates RNA from ADT/ATAC features internally; caller never needs to subset manually
3. `var_names_make_unique()` on the RNA subset only
4. Mitochondrial gene detection — auto-detects `MT-` (human) or `mt-` (mouse) prefix; falls back to `adata.var['gene_ids']`
5. Per-cell QC metrics via `sc.pp.calculate_qc_metrics()` on RNA only:
   - `n_genes_by_counts` — genes detected per cell
   - `total_counts` — total RNA UMI count per cell
   - `pct_counts_mt` — mitochondrial read percentage
6. Doublet detection via Scrublet on RNA only — adds `obs['doublet_score']` and `obs['predicted_doublet']`; fails gracefully
7. Cell filter application — removes cells failing any RNA-based threshold
8. MuData assembly — filtered cells applied to all modalities; QC obs columns live on `mdata["rna"]` only
9. Optional HTML report generation via `qc_report.py`

**Default thresholds**

| Parameter | Default | What it removes |
|-----------|---------|----------------|
| `min_genes` | 200 | Empty droplets / dead cells |
| `max_genes` | 6000 | Likely multiplets |
| `max_mt_pct` | 20% | Lysed / dying cells |
| `remove_doublets` | True | Scrublet-detected doublets |

All thresholds are configurable per-dataset.

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Raw counts in `adata.X`; may be mixed modality |
| `modality` | str | `"auto"` | `"auto"` \| `"rna"` \| `"cite"` \| `"multiome"` |
| `min_genes` | int | 200 | Minimum genes per cell |
| `max_genes` | int | 6000 | Maximum genes per cell |
| `max_mt_pct` | float | 20.0 | Maximum MT% per cell |
| `remove_doublets` | bool | True | Remove Scrublet doublets |
| `generate_report` | bool | False | Write HTML report |
| `report_path` | str | `"reports/qc_report.html"` | HTML output path |
| `sample_name` | str | `"sample"` | Label in report header |

**Input**

```
adata  : AnnData   →  raw counts in adata.X (output of ingest.py)
                       May contain mixed modalities (RNA + ADT or RNA + ATAC)
```

**Output**

```
mdata              : MuData   →  one AnnData per modality, filtered cells only
mdata["rna"]       : AnnData  →  RNA features, QC metrics in .obs (always present)
mdata["adt"]       : AnnData  →  ADT features, clean .obs (CITE-seq only)
mdata["atac"]      : AnnData  →  ATAC peaks, clean .obs (Multiome only)
metrics            : dict     →  cell counts, per-filter removal counts, medians,
                                 thresholds, detected modality
```

**Usage**

```python
from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc

# Works for plain RNA, CITE-seq, and Multiome — no manual subsetting needed
adata = load_dataset("data/benchmark/GSE194122_cite_raw_only.h5ad")

mdata, metrics = run_qc(adata)

# Access results
adata_rna = mdata["rna"]          # RNA AnnData with QC metrics in .obs
adata_adt = mdata["adt"]          # ADT AnnData (CITE-seq only)
print(metrics["modality"])        # "cite"
print(metrics["n_cells_output"])  # cells that passed QC

# Custom thresholds + HTML report
mdata, metrics = run_qc(
    adata,
    min_genes=300,
    max_genes=5000,
    max_mt_pct=15.0,
    generate_report=True,
    report_path="reports/qc_GSE194122.html",
    sample_name="GSE194122_BMMC",
)

# Pass RNA to normalization
from pipeline.modules.qc.normalize import normalize
adata_norm = normalize(mdata["rna"])
```

**Helper functions (importable)**

```python
from pipeline.modules.qc.qc import _detect_modality

modality = _detect_modality(adata)  # "rna" | "cite" | "multiome"
```

**Connects to**: `normalize.py` — pass `mdata["rna"]`

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
adata_raw      : AnnData  →  pre-filter RNA AnnData (after metrics computed)
adata_filtered : AnnData  →  post-filter RNA AnnData
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


---

## 5. `normalize.py`

**What it does**

Normalizes raw-count RNA AnnData (the `mdata["rna"]` slot from `qc.py`) and
selects highly variable genes for dimensionality reduction. Produces two layers
so every downstream module can access either raw or normalized counts without
recomputation.

**Why it exists**

Raw counts are not directly comparable across cells — cells with more RNA
captured will appear to express every gene more highly. Normalization removes
this cell-level sequencing depth effect. HVG selection then reduces the gene
space to the genes that carry the most biological signal, which speeds up PCA,
UMAP, and clustering without sacrificing accuracy.

**Steps performed (in order)**

1. Input validation — rejects non-AnnData inputs and already-normalized matrices
2. Save raw counts to `layers['counts']` before any modification
3. HVG selection on raw counts (`seurat_v3` flavor only — run pre-normalization
   as the method requires integer counts for its variance model)
4. Normalize per cell to `target_sum` counts (`sc.pp.normalize_total`)
5. log1p transform (`sc.pp.log1p`)
6. Save log1p-normalized values to `layers['logcounts']` (Seurat convention)
7. HVG selection for non-`seurat_v3` flavors (run post log1p)
8. Store all parameters and software versions in `uns['omicsage_normalization']`

**Layer layout after normalize()**

| Slot | Contents |
|------|----------|
| `.X` | log1p-normalized values (same as `logcounts`) |
| `layers['counts']` | Raw integer counts (preserved from input) |
| `layers['logcounts']` | log1p CP10K values — Seurat convention |
| `var['highly_variable']` | Boolean HVG flag (default: top 2000 genes) |
| `uns['omicsage_normalization']` | Full provenance record |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Raw counts in `.X` — pass `mdata["rna"]` from `run_qc()` |
| `target_sum` | float | `1e4` | Per-cell normalization target (CP10K) |
| `n_top_genes` | int | 2000 | Number of highly variable genes to select |
| `hvg_flavor` | str | `"seurat_v3"` | HVG method — `"seurat_v3"` \| `"seurat"` \| `"cell_ranger"` |
| `batch_key` | str | None | `obs` column for per-batch HVG selection (e.g. `"batch"`, `"DonorID"`) |
| `min_mean` | float | 0.0125 | HVG filter — min mean expression (non-`seurat_v3` only) |
| `max_mean` | float | 3.0 | HVG filter — max mean expression (non-`seurat_v3` only) |
| `min_disp` | float | 0.5 | HVG filter — min dispersion (non-`seurat_v3` only) |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**Note on `hvg_flavor`**: `seurat_v3` fits a regularized negative binomial
variance model on raw counts and requires `scikit-misc` (`pip install scikit-misc`).
It also requires ≥ ~200 cells to avoid numerical singularities — use `seurat`
flavor only for small test fixtures.

**Note on `batch_key`**: When provided, `sc.pp.highly_variable_genes` runs
per batch and flags a gene as highly variable if it qualifies in at least one
batch. Scanpy adds `var['highly_variable_nbatches']` to record how many batches
called each gene variable. Recommended for multi-donor / multi-site datasets.

**Input**

```
adata  : AnnData   →  raw counts in adata.X  (output of run_qc() → mdata["rna"])
```

**Output**

```
adata_norm.X                          →  log1p-normalized values
adata_norm.layers['counts']           →  raw integer counts (preserved)
adata_norm.layers['logcounts']        →  log1p CP10K values
adata_norm.var['highly_variable']     →  boolean HVG flag
adata_norm.var['highly_variable_nbatches']  →  per-batch HVG count (if batch_key used)
adata_norm.uns['omicsage_normalization']    →  provenance record
metrics  : dict                       →  n_cells, n_genes, n_hvg_selected,
                                          target_sum, hvg_flavor, batch_key,
                                          mean_counts_per_cell_after_norm,
                                          log1p_applied, raw_counts_in_layer,
                                          normalized_in_layer
```

**Usage**

```python
from pipeline.modules.qc.normalize import normalize

# Minimal — production defaults
adata_norm, metrics = normalize(mdata["rna"])

# With batch correction for multi-donor data
adata_norm, metrics = normalize(
    mdata["rna"],
    target_sum=1e4,
    n_top_genes=2000,
    hvg_flavor="seurat_v3",
    batch_key="batch",         # or "DonorID", "Site"
    inplace=False,
)

print(metrics["n_hvg_selected"])   # 2000
print(metrics["batch_key"])        # "batch"

# Access layers
adata_norm.layers["counts"]     # raw counts
adata_norm.layers["logcounts"]  # log1p normalized
```

**Connects to**: `normalization_report.py` for the report, then `reduce.py` for PCA + UMAP

---

## 6. `normalization_report.py`

**What it does**

Generates a self-contained HTML report for the normalization step.
Contains four figures and two summary tables. The HTML file is fully
portable — all figures are base64-embedded PNGs, no internet required.

**Why it exists**

Same philosophy as `qc_report.py`: every pipeline step produces a
human-readable output. The normalization report lets a biologist confirm
that HVG selection looks sensible (the scatter should show a clear
population of variable genes) and that library sizes are uniform after
normalization before dimensionality reduction begins.

**Report contents**

| Figure | What it shows |
|--------|--------------|
| HVG scatter | Mean expression vs normalised variance — HVGs in orange, background in grey |
| Library size distribution | Per-cell total counts before vs after normalization (violin) |
| Top 20 HVGs | Bar chart of the most variable genes by name |
| Gene detection rate | Cumulative fraction of cells each gene is detected in |

Plus a summary table (cells, genes, HVG count, flavor, batch key, target sum)
and a provenance table pulled directly from `uns['omicsage_normalization']`.

**Input**

```
adata_norm   : AnnData   →  output of normalize() — must have layers['counts']
                             and layers['logcounts']
metrics      : dict      →  output of normalize()
report_path  : str       →  where to write the .html file
dataset_name : str       →  label shown in the report header
```

**Output**

```
reports/output/normalization_report.html   →  self-contained HTML (~1-3 MB)
```

**Usage**

From notebook:

```python
from reports.normalization_report import run_normalization_report

run_normalization_report(
    adata_norm=adata_norm,
    metrics=metrics,
    report_path="reports/output/normalization_report.html",
    dataset_name="GSE194122_CITE",
)
```

From CLI (runs normalize + saves h5ad + generates report in one command):

```bash
python reports/normalization_report.py \
    --input  data/processed/GSE194122_cite_rna_qc.h5ad \
    --output data/processed/GSE194122_cite_normalized.h5ad \
    --report reports/output/normalization_report.html \
    --dataset GSE194122_CITE \
    --batch-key batch
```

**Dependencies**: `matplotlib` only — no Plotly, no Bokeh, no JS frameworks.

**Connects to**: `reduce.py` — pass `adata_norm` to PCA + UMAP


---

## Module Data Flow

```
Raw file (.h5ad / .h5 / MTX dir)
        │
        ▼
   ingest.py          → adata.X = raw counts (may contain mixed modalities)
        │
        ▼
     qc.py            → MuData: mdata["rna"] + mdata["adt"] / mdata["atac"]
        │                       │
        │                       ▼
        │               qc_report.py  → reports/qc_report.html
        │
        ├── mdata["rna"]
        │       │
        │       ▼
        │   normalize.py  → layers['counts'] + layers['logcounts'] + HVGs
        │       │                   │
        │       │                   ▼
        │       │       normalization_report.py → reports/output/normalization_report.html
        │       │
        │       ▼
        │   reduce.py     → obsm['X_pca'] + obsm['X_umap']   ← NEXT STEP
        │
        ├── mdata["adt"]  → ADT QC + CLR normalization (future phase)
        └── mdata["atac"] → ATAC QC (Phase 4)
```

---

## Running the Full QC Step

```bash
cd ~/OmicSage
conda activate omicsage

python3 -c "
from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc

# Pass the full mixed AnnData — run_qc() handles modality detection internally
adata = load_dataset('data/benchmark/GSE194122_cite_raw_only.h5ad')
mdata, metrics = run_qc(
    adata,
    generate_report=True,
    report_path='reports/qc_BMMC.html',
    sample_name='GSE194122_BMMC',
)
print('Modality detected:', metrics['modality'])
print('Cells before:     ', metrics['n_cells_input'])
print('Cells after:      ', metrics['n_cells_output'])
print('Removed:          ', metrics['n_cells_removed'])
print('MuData keys:      ', list(mdata.mod.keys()))
"
```

---

## Tests

| Test file | What it covers |
|-----------|---------------|
| `tests/test_phase0_structure.py` | Repo structure, imports, config schema |
| `tests/test_ingest.py` | Format detection, raw count extraction, all three loaders |
| `tests/test_qc.py` | MT detection, metric computation, filtering, Scrublet, ground-truth validation, modality detection, MuData structure, ADT/ATAC preservation |
| `tests/test_normalize.py` | Raw count preservation, normalization correctness, log1p, HVG selection, HVG count accuracy, batch_key, logcounts layer, provenance, mutation guard, input validation |

Run all QC + normalization tests:

```bash
conda activate omicsage
python -m pytest tests/test_phase0_structure.py tests/test_ingest.py tests/test_qc.py tests/test_normalize.py -v
# Expected: 42 passed (test_qc.py), 12 passed (test_normalize.py), 2 skipped
```
