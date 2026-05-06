# OmicSage — QC Module: Script Reference

> Location: `pipeline/modules/qc/`
> Phase: 1 — Core scRNA Pipeline
> Last updated: 2026-05-06

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

**Connects to**: `normalize.py` — pass `mdata["rna"]`

---

## 3. `qc_report.py`

**What it does**

Generates a self-contained HTML QC report after filtering. The report
is a single `.html` file with all plots embedded as base64 PNG images —
no external dependencies, no internet required, opens in any browser.

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

**Usage**

Called automatically by `run_qc()` when `generate_report=True`.

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

---

## 4. `data_report.py`

**What it does**

Generates a data intake HTML report summarising the contents of a raw
h5ad file *before any QC is applied*. Uses `backed='r'` mode — only
metadata is loaded, the full count matrix is never read. Works on files
of any size including large multiome files.

**Usage**

```bash
python pipeline/modules/qc/data_report.py \
  --input  data/benchmark/GSE194122_BMMC.h5ad \
  --geo    GSE194122 \
  --output reports/data_intake_BMMC.html
```

**Connects to**: Read the report → decide QC thresholds → run `qc.py`

---

## 5. `normalize.py`

**What it does**

Normalizes raw-count RNA AnnData (the `mdata["rna"]` slot from `qc.py`) and
selects highly variable genes for dimensionality reduction. Produces two layers
so every downstream module can access either raw or normalized counts without
recomputation.

**Steps performed (in order)**

1. Input validation — rejects non-AnnData inputs and already-normalized matrices
2. Save raw counts to `layers['counts']` before any modification
3. HVG selection on raw counts (`seurat_v3` flavor only — run pre-normalization)
4. Normalize per cell to `target_sum` counts (`sc.pp.normalize_total`)
5. log1p transform (`sc.pp.log1p`)
6. Save log1p-normalized values to `layers['logcounts']` (Seurat convention)
7. HVG selection for non-`seurat_v3` flavors (run post log1p)
8. Store all parameters, software versions, and timestamp in `uns['omicsage_normalization']`

**Layer layout after normalize()**

| Slot | Contents |
|------|----------|
| `.X` | log1p-normalized values (same as `logcounts`) |
| `layers['counts']` | Raw integer counts (preserved from input) |
| `layers['logcounts']` | log1p CP10K values — Seurat convention |
| `var['highly_variable']` | Boolean HVG flag (default: top 2000 genes) |
| `uns['omicsage_normalization']` | Full provenance record including timestamp |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Raw counts in `.X` — pass `mdata["rna"]` from `run_qc()` |
| `target_sum` | float | `1e4` | Per-cell normalization target (CP10K) |
| `n_top_genes` | int | 2000 | Number of highly variable genes to select |
| `hvg_flavor` | str | `"seurat_v3"` | HVG method — `"seurat_v3"` \| `"seurat"` \| `"cell_ranger"` |
| `batch_key` | str | None | `obs` column for per-batch HVG selection |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**Usage**

```python
from pipeline.modules.qc.normalize import normalize

adata_norm, metrics = normalize(
    mdata["rna"],
    target_sum=1e4,
    n_top_genes=2000,
    hvg_flavor="seurat_v3",
    batch_key="batch",
    inplace=False,
)
```

**Connects to**: `normalization_report.py` for the report, then `reduce.py` for PCA + UMAP

---

## 6. `normalization_report.py`

**What it does**

Generates a self-contained HTML report for the normalization step.
Contains four figures and two summary tables. Callable from notebook or CLI.

**Report contents**

| Figure | What it shows |
|--------|--------------|
| HVG scatter | Mean expression vs normalised variance — HVGs in orange, background in grey |
| Library size distribution | Per-cell total counts before vs after normalization (violin) |
| Top 20 HVGs | Bar chart of the most variable genes by name |
| Gene detection rate | Cumulative fraction of cells each gene is detected in |

**Usage**

```python
from reports.normalization_report import run_normalization_report

run_normalization_report(
    adata_norm=adata_norm,
    metrics=metrics,
    report_path="reports/output/normalization_report.html",
    dataset_name="GSE194122_CITE",
)
```

```bash
python reports/normalization_report.py \
    --input  data/processed/GSE194122_cite_rna_qc.h5ad \
    --output data/processed/GSE194122_cite_normalized.h5ad \
    --report reports/output/normalization_report.html \
    --dataset GSE194122_CITE \
    --batch-key batch
```

**Connects to**: `reduce.py` — pass `adata_norm` to PCA + UMAP

---

## 7. `reduce.py`

**What it does**

Runs PCA, neighbor graph construction, UMAP (and optionally t-SNE) on a
normalized AnnData (output of `normalize.py`). The number of PCs used for
the neighbor graph is chosen automatically by default using elbow detection
on the variance-explained curve (`kneed`). The user can override with
`n_pcs=<int>`.

**Why it exists**

Raw gene expression space has thousands of dimensions — most of which are
noise. PCA compresses this into a low-dimensional space capturing the major
axes of biological variation. UMAP then projects this into 2D for visualization.
The neighbor graph (kNN) built on the PCA embedding is the foundation for
Leiden clustering in the next step.

Choosing n_pcs data-adaptively (rather than hardcoding 30) ensures the neighbor
graph captures real biological variation without including noise PCs that would
blur cluster boundaries.

**Steps performed (in order)**

1. Input validation — warns if raw counts are passed instead of normalized data
2. PCA on HVG subset only (`use_highly_variable=True`, `svd_solver='arpack'`)
3. Auto-select n_pcs via elbow detection (kneed) — fallback chain:
   elbow → variance threshold (85%) → fixed (30) — each step logs a warning
4. Neighbor graph (`sc.pp.neighbors`, `n_neighbors=15`, `n_pcs=selected`)
5. UMAP (`sc.tl.umap`) — always computed
6. t-SNE (`sc.tl.tsne`) — optional, default off (slow on large datasets)
7. Store all parameters, software versions, and timestamp in `uns['omicsage_reduce']`

**obsm / obsp layout after reduce()**

| Slot | Contents |
|------|----------|
| `obsm['X_pca']` | PCA embedding (n_cells × n_comps) |
| `obsm['X_umap']` | UMAP embedding (n_cells × 2) |
| `obsm['X_tsne']` | t-SNE embedding (n_cells × 2) — only if `run_tsne=True` |
| `obsp['connectivities']` | Sparse kNN connectivity matrix |
| `obsp['distances']` | Sparse kNN distance matrix |
| `uns['omicsage_reduce']` | Full provenance record including timestamp |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Normalized AnnData — log1p in `.X`, HVGs in `var['highly_variable']` |
| `n_comps` | int | 50 | Number of PCA components to compute |
| `n_pcs` | int | None | PCs to use for neighbor graph — None = auto-select |
| `n_pcs_method` | str | `"elbow"` | `"elbow"` \| `"variance"` \| `"fixed"` |
| `variance_threshold` | float | 0.85 | Cumulative variance target for `n_pcs_method="variance"` |
| `n_neighbors` | int | 15 | Number of neighbors for kNN graph |
| `run_tsne` | bool | False | Also compute t-SNE (slow on large datasets) |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |
| `random_state` | int | 0 | Reproducibility seed for all stochastic steps |

**PC selection methods**

| Method | How it works | When to use |
|--------|-------------|-------------|
| `"elbow"` (default) | kneed KneeLocator on variance-explained curve; falls back to `"variance"` | Most datasets |
| `"variance"` | Keep PCs until cumulative variance ≥ `variance_threshold` | Flat scree curves |
| `"fixed"` | Always use `min(30, n_comps)` | Reproducibility / debugging |
| manual `n_pcs=N` | Bypasses auto-selection entirely | Expert override |

**Input**

```
adata  : AnnData   →  normalized AnnData (output of normalize())
                       Must have log1p values in .X and HVGs flagged in var['highly_variable']
```

**Output**

```
adata_reduced.obsm['X_pca']     →  PCA embedding
adata_reduced.obsm['X_umap']    →  UMAP embedding
adata_reduced.obsm['X_tsne']    →  t-SNE embedding (if run_tsne=True)
adata_reduced.obsp['connectivities']  →  kNN graph
adata_reduced.obsp['distances']       →  kNN distances
adata_reduced.uns['omicsage_reduce']  →  provenance record
metrics  : dict                 →  n_cells, n_genes, n_hvg_used, n_comps_computed,
                                    n_pcs_used, pc_selection_method,
                                    cumulative_variance_explained,
                                    variance_explained_per_pc, n_neighbors,
                                    run_tsne, embeddings_computed
```

**Usage**

```python
from pipeline.modules.qc.reduce import reduce

# Production defaults — elbow PC selection
adata_reduced, metrics = reduce(adata_norm)

# Manual PC override
adata_reduced, metrics = reduce(adata_norm, n_pcs=20)

# With t-SNE and variance-based PC selection
adata_reduced, metrics = reduce(
    adata_norm,
    n_comps=50,
    n_pcs=None,
    n_pcs_method="variance",
    variance_threshold=0.85,
    n_neighbors=15,
    run_tsne=True,
    inplace=False,
)

print(metrics["n_pcs_used"])            # e.g. 18
print(metrics["pc_selection_method"])   # "elbow"
print(metrics["cumulative_variance_explained"])  # e.g. 0.81
```

**Connects to**: `reduce_report.py` for the report, then `cluster.py` for Leiden clustering

---

## 8. `reduce_report.py`

**What it does**

Generates a self-contained HTML report for the dimensionality reduction step.
Contains five figures and two summary tables. Callable from notebook or CLI.

**Why it exists**

The scree plot is the primary tool for validating PC selection. The QC overlays
on UMAP detect whether major axes of variation are driven by biology or technical
factors (sequencing depth, MT%, doublets). If UMAP structure is dominated by
`total_counts` or `pct_counts_mt`, covariates should be regressed out before
clustering.

**Report contents**

| Figure | What it shows |
|--------|--------------|
| Scree plot (2 panels) | Per-PC variance bars + cumulative curve; red dashed line at selected n_pcs |
| UMAP × QC (2×2 grid) | n_genes, total_counts, MT%, doublet_score overlaid on UMAP |
| PCA × QC (1×2) | PC1 vs PC2 coloured by n_genes and total_counts |
| UMAP × batch (optional) | One colour per batch — only rendered when `batch_key` is passed |

**Usage**

```python
from reports.reduce_report import run_reduce_report

run_reduce_report(
    adata_reduced=adata_reduced,
    metrics=metrics,
    report_path="reports/output/reduce_report.html",
    dataset_name="GSE194122_CITE",
    batch_key="batch",    # optional — adds batch UMAP panel
)
```

```bash
python reports/reduce_report.py \
    --input  data/processed/GSE194122_cite_normalized.h5ad \
    --output data/processed/GSE194122_cite_reduced.h5ad \
    --report reports/output/reduce_report.html \
    --dataset GSE194122_CITE \
    --batch-key batch
```

**Connects to**: `cluster.py` — pass `adata_reduced` to Leiden clustering

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
        │   reduce.py     → obsm['X_pca'] + obsm['X_umap'] + obsp['connectivities']
        │       │                   │
        │       │                   ▼
        │       │       reduce_report.py → reports/output/reduce_report.html
        │       │
        │       ▼
        │   cluster.py    → obs['leiden_*']    ← NEXT STEP
        │
        ├── mdata["adt"]  → ADT QC + CLR normalization (future phase)
        └── mdata["atac"] → ATAC QC (Phase 4)
```

---

## Tests

| Test file | What it covers |
|-----------|---------------|
| `tests/test_phase0_structure.py` | Repo structure, imports, config schema |
| `tests/test_ingest.py` | Format detection, raw count extraction, all three loaders |
| `tests/test_qc.py` | MT detection, metric computation, filtering, Scrublet, ground-truth validation, modality detection, MuData structure, ADT/ATAC preservation |
| `tests/test_normalize.py` | Raw count preservation, normalization correctness, log1p, HVG selection, HVG count accuracy, batch_key, logcounts layer, provenance, mutation guard, input validation |
| `tests/test_reduce.py` | PCA/UMAP shapes, HVG-only PCA, neighbor graph, provenance, inplace guard, t-SNE optional, elbow/variance/manual PC selection, no-HVG fallback |

Run all tests:

```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: 162 passed (approximate), 2 skipped
```

Run just the reduce tests:

```bash
python -m pytest tests/test_reduce.py -v
# Expected: 12 passed
```
