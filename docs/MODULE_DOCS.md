# OmicSage — Module Reference

> Locations: `pipeline/modules/qc/` · `pipeline/modules/annotation/` · `reports/`
> Phase: 2 — Report Engine
> Last updated: 2026-05-12 (session 8)

This document describes every script in the pipeline — what it does,
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

**Steps performed (in order)**

1. Input validation — warns if raw counts are passed instead of normalized data
2. PCA on HVG subset only (`use_highly_variable=True`, `svd_solver='arpack'`)
3. Auto-select n_pcs via elbow detection (kneed) — fallback chain: elbow → variance threshold (85%) → fixed (30)
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
```

**Connects to**: `reduce_report.py` for the report, then `cluster.py` for Leiden clustering

---

## 8. `reduce_report.py`

**What it does**

Generates a self-contained HTML report for the dimensionality reduction step.
Contains five figures and two summary tables. Callable from notebook or CLI.

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

**Connects to**: `cluster.py` — pass `adata_reduced` to Leiden clustering

---

## 9. `cluster.py`

**What it does**

Runs Leiden community detection at multiple resolutions on the kNN graph
produced by `reduce.py`. Computes a silhouette score per resolution to
quantify cluster separation, and either auto-selects the best resolution or
uses a caller-supplied override.

**Why it exists**

Choosing a clustering resolution is one of the most consequential decisions
in single-cell analysis. Too low → biologically distinct populations are merged.
Too high → one real population gets split arbitrarily. Sweeping multiple
resolutions and scoring each with silhouette makes the trade-off explicit and
reproducible. The `best_resolution_override` parameter lets the analyst
over-cluster deliberately before annotation — the standard approach is to
over-cluster slightly then merge during annotation, rather than to split.

**Steps performed (in order)**

1. Input validation — checks `obsm['X_pca']` and `obsp['connectivities']` are present
2. Resolution sweep — runs `sc.tl.leiden` at each resolution in `resolution_range`
3. Stores each result in `obs[f'leiden_{res}']`
4. Computes silhouette score per resolution on X_pca (subsampled above 10k cells)
5. Selects best resolution — override if supplied, otherwise argmax of silhouette scores
6. Copies best-resolution labels to `obs['leiden']` as a convenience key
7. Stores provenance in `uns['omicsage_cluster']` with string-keyed dicts (HDF5 requirement)

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Reduced AnnData — must have `obsm['X_pca']` and `obsp['connectivities']` |
| `resolution_range` | list[float] | `[0.2, 0.4, 0.6, 0.8, 1.0]` | Leiden resolutions to sweep |
| `best_resolution_override` | float | None | Pin a specific resolution; None = silhouette auto-select |
| `pca_key` | str | `"X_pca"` | obsm key for silhouette scoring |
| `connectivities_key` | str | `"connectivities"` | obsp key for Leiden adjacency |
| `random_state` | int | 0 | Reproducibility seed |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**Output**

```
adata.obs[f'leiden_{res}']    →  cluster labels at each resolution tested
adata.obs['leiden']           →  labels at the selected resolution (convenience key)
adata.uns['omicsage_cluster'] →  provenance record
metrics : dict                →  resolutions, n_clusters, silhouette_scores,
                                  best_resolution, best_silhouette, best_n_clusters
```

**Usage**

```python
from pipeline.modules.qc.cluster import cluster

# Silhouette auto-select
adata_clustered, metrics = cluster(adata_reduced)

# Override resolution (recommended before annotation)
adata_clustered, metrics = cluster(
    adata_reduced,
    resolution_range=[0.2, 0.4, 0.6, 0.8, 1.0],
    best_resolution_override=0.8,
)
```

**Important note on resolution selection**

Silhouette score optimises for geometric separation in PCA space and tends to
select low resolutions. For annotation purposes it is usually better to
over-cluster and merge than to under-cluster and split. Use
`best_resolution_override` to pin a higher resolution when the dataset has many
expected cell types (e.g. PBMC / BMMC datasets).

**Connects to**: `cluster_report.py` for the report, then `annotate.py` for cell type annotation

---

## 10. `cluster_report.py`

**What it does**

Generates a self-contained HTML report for the Leiden clustering step.

**Report contents**

| Figure | What it shows |
|--------|--------------|
| UMAP grid (one panel per resolution) | Cluster assignments per resolution; gold border on the selected panel |
| Silhouette bar chart | Score per resolution; orange bar = selected; value labels on each bar |
| Cluster size distribution | Cell counts per cluster at best resolution, sorted descending; median line |
| UMAP × ground-truth cell_type (optional) | Only rendered when `obs['cell_type']` is present |

**Usage**

```python
from reports.cluster_report import run_cluster_report

run_cluster_report(
    adata_clustered=adata_clustered,
    metrics=metrics,
    report_path="reports/output/cluster_report.html",
    dataset_name="GSE194122_CITE",
)
```

**Connects to**: `annotate.py` — pass `adata_clustered`

---

## 11. `annotate.py`
 
**What it does**
 
Assigns cell type labels to Leiden clusters using up to five automated
methods, then combines them into a consensus label via weighted majority
vote with semantic label normalisation.
All methods operate at the cluster level (not per cell), so every cell
in the same Leiden cluster receives the same consensus label.
 
**Why it exists**
 
No single annotation method is universally reliable. CellTypist is strong
on immune cell types but can conflate closely related states. Marker gene
scoring is interpretable but depends on the quality of the gene list.
ScType uses a curated tissue-specific database. SingleR correlates against
a bulk RNA-seq pseudobulk reference. Running multiple methods and taking a
weighted majority vote produces more robust annotations than any one method
alone — and makes method disagreements explicit through the confidence score.
 
---
 
### Methods
 
| Method | Key | How it works | Reference column(s) |
|--------|-----|-------------|---------------------|
| CellTypist | `"celltypist"` | Per-cell prediction (any model), then cluster-level majority vote done by OmicSage. Each model produces one `celltypist_<stem>` column. `Immune_All_High.pkl` → `celltypist_coarse`; `Immune_All_Low.pkl` → `celltypist_fine`; any other model → `celltypist_<model_stem>`. | `celltypist_*` (one per model) |
| Marker scoring | `"markers"` | Mean log-normalised expression of each marker gene set per cluster; assigns the highest-scoring type. | `cell_type_markers` |
| ScType | `"sctype"` | Positive/negative marker scoring against ScTypeDB; fetched fresh from GitHub at runtime unless `sctype_db_path` is set. | `cell_type_sctype` |
| SingleR | `"singler"` | Spearman correlation against a pseudobulk reference via the `singler` Python package (C++ bindings). Per-cell labels; delta score used for pruning. | `cell_type_singler`, `singler_delta` |
| scANVI | `"scanvi"` | Transfer learning from a pre-trained scANVI model. Posterior probability used as fractional vote weight. Skipped in CI (`OMICSAGE_CI=true`). | `cell_type_scanvi` |
| Majority vote | `"vote"` | Weighted consensus across all methods that ran. Requires at least `"celltypist"` or `"markers"` to also be in methods. | `cell_type_vote`, `cell_type_confidence` |
 
---
 
### Steps performed (in order)
 
1. Validate methods and prerequisites (`"vote"` requires `"celltypist"` or `"markers"`)
2. Preserve ground-truth: copy `obs['cell_type']` → `obs['cell_type_groundtruth']` before any writes
3. Run CellTypist — builds a CP10K log1p copy, downloads requested models to `celltypist_models_dir`, runs per-cell prediction then applies OmicSage's own cluster-level majority vote; writes one `celltypist_*` column per model
4. Run marker scoring — mean expression per cluster per marker set; writes `obs['cell_type_markers']`
5. Run ScType — fetches ScTypeDB from GitHub (or loads from `sctype_db_path`); writes `obs['cell_type_sctype']`
6. Run SingleR — loads the requested celldex reference (cached on first use); writes `obs['cell_type_singler']` and `obs['singler_delta']`
7. Run scANVI — loads a pre-trained model; writes `obs['cell_type_scanvi']`
8. Run majority vote — collects all active method columns, normalises labels, accumulates weighted votes, elects winner; writes `obs['cell_type_vote']` and `obs['cell_type_confidence']`
9. Write provenance to `uns['omicsage_annotate']`
Each step is wrapped in try/except — if a method fails it warns and continues, so a network error or missing package never crashes the pipeline.
 
---
 
### How the majority vote works
 
**Stage 1 — collect labels.** For each cluster and each active method, take the majority label among all cells in that cluster.
 
**Stage 2 — normalise.** Each raw label is passed through `_normalise_label()` to produce a canonical root concept. This ensures synonymous labels from different methods group together as the same vote rather than splitting:
 
```
"Late erythroid"                         →  "erythroid"
"Erythroid cells"                        →  "erythroid"
"Erythroid-like and erythroid precursor" →  "erythroid"
"Classical monocytes"                    →  "monocyte"
"Monocytes"                              →  "monocyte"
"CD4+ T cells"                           →  "t"
"T_cell"                                 →  "t"
```
 
Strips applied: underscores → parentheses → stage prefixes (`Late`, `Classical`, `Naive`, `CD4+`, …) → `pre-`/`pro-` → suffixes (`-like and …`, `precursor`, `cells`, `cell`) → singularise trailing `s`.
 
**Stage 3 — accumulate weighted votes.** Each method contributes its weight to the canonical group its label maps to:
 
| Source | Weight |
|--------|--------|
| Each `celltypist_*` column | 1 |
| `cell_type_markers` | 1 |
| `cell_type_sctype` | 1 (+1 bonus for parenchymal types: Hepatocyte, Fibroblast, Endothelial) |
| `cell_type_singler` | 1 |
| `cell_type_scanvi` | `mean_posterior × 2` (fractional, 0–2) |
 
**Stage 4 — elect winner and compute confidence.**
The canonical group with the highest accumulated weight wins.
`confidence = winning_weight / total_weight`.
All four integer-weight methods agreeing on one root concept → `confidence = 1.0`, even if the raw strings differ.
 
**Tiebreak** (when two canonical groups share the highest weight):
priority order is `celltypist_fine → singler → markers → sctype → celltypist_coarse → other celltypist_*`.
 
**Stage 5 — resolve reported label.**
The final `cell_type_vote` is the **raw label from the highest-priority source** within the winning canonical group — not the lowercase canonical form — to preserve the most informative string (e.g. `"Late erythroid"` not `"erythroid"`).
 
---
 
### obs columns written
 
| Column | Source | Notes |
|--------|--------|-------|
| `celltypist_<stem>` | CellTypist, one per model | e.g. `celltypist_coarse`, `celltypist_fine`, `celltypist_Healthy_COVID19_PBMC` |
| `cell_type_markers` | Marker gene scoring | Best-scoring type per cluster |
| `cell_type_sctype` | ScType | Best ScType label per cluster |
| `cell_type_singler` | SingleR | Per-cell label; `"Unassigned"` for low-delta cells |
| `singler_delta` | SingleR | Per-cell delta score (top − second correlation); NaN for `"Unassigned"` |
| `cell_type_scanvi` | scANVI | Per-cell transfer label |
| `cell_type_groundtruth` | Copy of original `obs['cell_type']` | Only written if column exists before annotation |
| `cell_type_vote` | Majority vote | Primary consensus label — use this in all downstream steps |
| `cell_type_confidence` | Majority vote | 0.0–1.0; flag clusters below 0.5 for manual review |
 
---
 
### Built-in marker sets (`MARKER_SETS`)
 
Covers immune (BMMC) and parenchymal (HCC) cell types:
 
| Category | Cell types |
|----------|-----------|
| Myeloid | Macrophage, Monocyte, DC, pDC, Mast_cell |
| Lymphoid | T_cell, CD8_T_cell, NK_ILC, B_cell, Plasma_cell |
| Progenitors / erythroid | Progenitor, Erythroid, Platelet |
| Non-immune / parenchymal | Hepatocyte, Fibroblast, Endothelial |
 
Pass a custom dict via `marker_sets={"MyType": ["GENE1", "GENE2"]}` to override or extend.
 
---
 
### Parameters
 
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Clustered AnnData — must have `obs[leiden_col]` and `layers['logcounts']` |
| `methods` | list[str] | `["celltypist", "markers", "sctype", "singler", "vote"]` | Methods to run |
| `leiden_col` | str | `"leiden"` | obs column with cluster labels |
| `celltypist_models` | list[str] | `["Immune_All_High.pkl", "Immune_All_Low.pkl"]` | Any CellTypist model filenames; each produces one `celltypist_*` obs column |
| `celltypist_models_dir` | Path | `data/references/celltypist/` | Local model cache directory |
| `marker_sets` | dict | `MARKER_SETS` | `{cell_type: [gene_list]}` |
| `tissue` | str | `"Immune system"` | ScType tissue filter. Valid: Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus |
| `sctype_db_path` | str / Path | `None` | Local ScTypeDB `.xlsx` file. `None` → fetch from GitHub each run |
| `scanvi_model` | str | `None` | Path to pre-trained scANVI model directory |
| `singler_ref` | str / Path | `None` | celldex reference name or file path. `None` → HPCA. Named options: `"hpca"`, `"blueprint_encode"`, `"dice"`, `"monaco_immune"`, `"novershtern_hematopoietic"`, `"mouse_rnaseq"`, `"hca"` |
| `singler_ref_label_col` | str | `"cell_type"` | obs column for cell type labels in a user-supplied H5AD reference |
| `inplace` | bool | `False` | Modify input AnnData in place; default makes a copy |
 
---
 
### Input / Output
 
```
Input
-----
adata : AnnData   →  clustered AnnData (output of cluster())
                     Must have obs['leiden'], layers['logcounts']
 
Output
------
adata_ann.obs['celltypist_*']          →  one column per CellTypist model run
adata_ann.obs['cell_type_markers']     →  marker-score label per cell
adata_ann.obs['cell_type_sctype']      →  ScType label per cell
adata_ann.obs['cell_type_singler']     →  SingleR label per cell
adata_ann.obs['singler_delta']         →  SingleR delta score per cell
adata_ann.obs['cell_type_scanvi']      →  scANVI transfer label per cell
adata_ann.obs['cell_type_groundtruth'] →  preserved ground-truth (if existed)
adata_ann.obs['cell_type_vote']        →  consensus label per cell  ← use downstream
adata_ann.obs['cell_type_confidence']  →  confidence score per cell (0.0–1.0)
adata_ann.uns['omicsage_annotate']     →  full provenance record
ann_dict : dict                        →  marker_score_df, vote_df, provenance
```
 
---
 
### Usage
 
```python
from pipeline.modules.annotation.annotate import annotate, MARKER_SETS
 
# Default — all five methods, HPCA reference for SingleR
adata_annotated, ann_dict = annotate(adata_clustered)
 
# BMMC-optimised — bone marrow SingleR reference, third CellTypist model
adata_annotated, ann_dict = annotate(
    adata_clustered,
    methods=["celltypist", "markers", "sctype", "singler", "vote"],
    celltypist_models=[
        "Immune_All_High.pkl",
        "Immune_All_Low.pkl",
        "Healthy_COVID19_PBMC.pkl",
    ],
    singler_ref="novershtern_hematopoietic",
)
 
# Offline run — pin ScTypeDB locally, skip SingleR network call
adata_annotated, ann_dict = annotate(
    adata_clustered,
    methods=["celltypist", "markers", "sctype", "vote"],
    sctype_db_path="data/references/ScTypeDB_full.xlsx",
)
 
# Custom marker sets (e.g. liver dataset)
liver_markers = {
    "Hepatocyte":  ["ALB", "APOC3", "TTR"],
    "T_cell":      ["CD3D", "CD3E"],
    "Endothelial": ["PECAM1", "VWF"],
}
adata_annotated, ann_dict = annotate(
    adata_clustered,
    methods=["markers"],
    marker_sets=liver_markers,
)
 
# Inspect consensus
print(adata_annotated.obs["cell_type_vote"].value_counts())
print(adata_annotated.obs["cell_type_confidence"].describe())
print(ann_dict["vote_df"][["n_cells", "final_label", "confidence"]])
```
 
---
 
### Important implementation notes
 
- **CellTypist now uses `majority_voting=False`** — OmicSage gets per-cell predictions and does its own cluster-level majority vote. This makes all CellTypist models behave identically; there is no dependency on CellTypist's internal `over_clustering` logic.
- Each CellTypist model produces its own `celltypist_<stem>` obs column and feeds independently into the vote. Running three models gives three columns and three votes.
- **`celltypist_coarse` and `celltypist_fine` only exist if those specific models were run.** No placeholder columns are written for models that weren't requested.
- celldex references for SingleR are downloaded once and cached at `~/.cache/omicsage/singler/<ref>.h5ad`. The first call per reference takes a few seconds; subsequent calls load from disk.
- ScType fetches `ScTypeDB_full.xlsx` from GitHub on every run when `sctype_db_path=None`. Set `sctype_db_path` after the first run to pin a version and enable offline use.
- The `CELLTYPIST_FOLDER` environment variable is set temporarily and then restored so the caller's environment is not polluted.
- Every method is wrapped in try/except — a failure warns and skips, so a network error or missing package never crashes the pipeline.
- `cell_type_vote` (not `cell_type`) is the output column — this avoids overwriting any existing ground-truth `cell_type` column from the dataset.
**Connects to**: `annotate_report.py` for the HTML report, then `deg.py` for differential expression

---

## 12. `annotate_report.py`

**What it does**

Generates a self-contained HTML report for the cell type annotation step.
Contains UMAP visualisations, per-cluster summary tables, and confidence
distribution plots.

**Report contents**

| Section | What it shows |
|---------|--------------|
| Summary cards | Methods run, clusters annotated, unique cell types, median confidence |
| UMAP × consensus vote | Cells coloured by `cell_type_vote` |
| UMAP × CellTypist fine | Cells coloured by `celltypist_fine` (for comparison) |
| Confidence distribution | Histogram of `cell_type_confidence` across all cells |
| Per-cluster table | Cluster → final label, confidence, n_cells, all method labels |

**Usage**

```python
from pipeline.modules.annotation.annotate_report import run_annotate_report

run_annotate_report(
    adata_annotated=adata_annotated,
    annotation_dict=ann_dict,
    report_path="reports/annotate_report.html",
    dataset_name="GSE194122_CITE",
)
```

**Connects to**: `deg.py` — pass `adata_annotated`

---

## 13. `deg.py`

**What it does**

Runs Wilcoxon rank-sum differential expression analysis across all annotated
cell types, one-vs-rest. Returns both a tidy per-group DataFrame and a
summary of the top 5 genes per group. Also supports optional pairwise
comparisons between specific groups.

**Why it exists**

Identifying genes that distinguish each cell type from all others is the
foundation for validating annotations (canonical marker genes should appear),
understanding cell biology, and — in the next step — pathway enrichment.
Performing DEG at the annotated cell-type level (rather than at the Leiden
cluster level) produces biologically interpretable results that connect
directly to the literature.

**Steps performed (in order)**

1. Resolve groupby column — tries `obs['cell_type_vote']`, falls back to `obs['leiden']` with a `UserWarning`
2. Warn about groups with < 10 cells (unreliable Wilcoxon statistics)
3. Set `adata.X = layers['logcounts']` when `use_raw=False` (default)
4. Run `sc.tl.rank_genes_groups` — Wilcoxon, `rankby_abs=True`, BH correction, `pts=True`
5. Extract results into tidy DataFrames (one per group): `gene, score, pval, logfc, pval_adj`
6. Apply logFC and pval_adj thresholds to filter to significant DEGs
7. Apply `exclude_gene_prefixes` post-filter (RPL/RPS/MT- etc.) — excluded genes still used in fold-change computation
8. Optionally run pairwise comparisons for specified group pairs
9. Build `summary_df` — top 5 significant DEGs per group in long format
10. Write provenance to `uns['omicsage_deg']`

**Why `rankby_abs=True` matters**

Without this flag, `sc.tl.rank_genes_groups` returns only the top-ranked
upregulated genes per group (ranked by raw score, always positive). Setting
`rankby_abs=True` ranks by |score|, so both strongly upregulated and strongly
downregulated genes appear in the results. This is required for volcano plots
to show both directions and for GSEA `direction="down"` or `direction="both"`
to have any genes to work with.

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Annotated AnnData — must have `layers['logcounts']` and `obs[groupby]` |
| `groupby` | str | `"cell_type_vote"` | obs column to group cells by |
| `leiden_col` | str | `"leiden"` | Fallback obs column if `groupby` is not present |
| `method` | str | `"wilcoxon"` | DEG method — `"wilcoxon"` \| `"t-test"` \| `"logreg"` |
| `min_logfc` | float | `0.25` | Minimum absolute log₂ fold-change threshold |
| `max_pval_adj` | float | `0.05` | Maximum BH-corrected adjusted p-value threshold |
| `n_genes` | int | `500` | Top genes to compute per group (before threshold filtering) |
| `exclude_gene_prefixes` | list[str] | None | Gene prefixes to remove from results e.g. `["RPL", "RPS", "MT-"]` |
| `use_raw` | bool | False | Use `adata.raw` / `adata.X` instead of `layers['logcounts']` |
| `pairwise_groups` | list[tuple] | None | List of `(group_a, group_b)` tuples for pairwise DEG |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**Output**

```
adata_deg.uns['omicsage_deg']           →  provenance record (includes exclude_gene_prefixes, rankby_abs)
adata_deg.uns['rank_genes_groups']      →  full scanpy output (all genes, all groups)

deg_dict['results']    : dict           →  {group: DataFrame(gene, score, pval, logfc, pval_adj)}
                                            — filtered to significant DEGs only; prefix-excluded genes removed
deg_dict['summary_df'] : DataFrame      →  long-format top 5 DEGs per group
                                            columns: group, rank, gene, logfc, pval_adj
deg_dict['provenance'] : dict           →  same as uns['omicsage_deg']
deg_dict['pairwise']   : dict           →  {(a, b): DataFrame} — only if pairwise_groups supplied
```

**Usage**

```python
from pipeline.modules.downstream.deg import deg

# Standard run — both directions, ribosomal genes filtered
adata_deg, deg_dict = deg(
    adata_annotated,
    groupby="cell_type_vote",
    method="wilcoxon",
    min_logfc=0.25,
    max_pval_adj=0.05,
    n_genes=500,
    exclude_gene_prefixes=["RPL", "RPS", "MT-"],
    inplace=False,
)

# Inspect significant DEGs for T cells
print(deg_dict["results"]["T_cell"].head(10))

# Top 5 DEGs per group (summary table)
print(deg_dict["summary_df"])

# Pairwise comparison
adata_deg, deg_dict = deg(
    adata_annotated,
    pairwise_groups=[("T_cell", "B_cell"), ("Monocyte", "DC")],
)
print(deg_dict["pairwise"][("T_cell", "B_cell")].head())
```

**Important implementation notes**

- `rankby_abs=True` is hardcoded — this is a deliberate design decision. Without it the volcano plots only show upregulated genes, which is misleading and breaks GSEA `direction="down"`.
- Threshold filtering happens *after* extraction, not inside `rank_genes_groups`. This means `uns['rank_genes_groups']` always contains the full `n_genes` ranked list, while `deg_dict['results']` contains only significant hits.
- `exclude_gene_prefixes` is applied *after* thresholding. Excluded genes still participate in the `rank_genes_groups` computation — removing them before would bias fold-changes for remaining genes.
- `pts=True` is always passed — adds fraction-expressing data to `uns['rank_genes_groups']` at no cost, used by dot plots.
- `n_genes=500` default prevents artificial capping. With the old 200 default, well-separated cell types (BMMC CITE-seq) returned exactly 200 genes because all computed genes passed thresholds.

**Connects to**: `deg_report.py` for the report, then `gsea.py` for pathway enrichment

---

## 14. `deg_report.py`

**What it does**

Generates a self-contained HTML report from `deg()` output. All figures
are embedded as base64 PNGs — no external files, opens in any browser.
Sections are built independently so one failed plot never kills the report.

**Report contents**

| Section | What it shows |
|---------|--------------|
| Run summary | Stat cards (groups, method, thresholds, n_genes computed, total significant DEGs) + per-group DEG counts table + prefix exclusion note (when applicable) |
| Top DEGs per group | Rowspan table — top 5 significant genes per group, ranked by adj. p-value; Direction column (▲ Up / ▼ Down) + log₂FC coloured red/blue |
| Volcano plots | One plot per group (default max 20, sorted by DEG count); up/down/NS classification; top N genes labelled; threshold dashed lines; visible note when truncation occurs |
| Dot plot | `sc.pl.dotplot` — dot size = fraction expressing, colour = mean expression — top 5 DEGs per group |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | AnnData returned by `deg()` |
| `deg_dict` | dict | required | `deg_dict` returned by `deg()` |
| `output_path` | str | `"reports/output/deg_report.html"` | HTML output path |
| `top_n_volcano` | int | 10 | Genes to label on each volcano plot |
| `top_n_dotplot` | int | 5 | Top DEGs per group to include in dot plot |
| `max_volcano_groups` | int | 20 | Max volcano plots rendered; when exceeded, groups sorted by DEG count and a visible note lists excluded groups |

**Usage**

```python
from reports.deg_report import generate_deg_report

report_path = generate_deg_report(
    adata=adata_deg,
    deg_dict=deg_dict,
    output_path="reports/deg_report.html",
)
print(f"Report → {report_path}")
```

**Connects to**: `gsea.py` — pass `deg_dict['results']` to pathway enrichment

---

## 15. `gsea.py`

**What it does**

Runs over-representation analysis (ORA) on DEG results using `gseapy.enrichr`.
For each cell type group, builds a directional query gene list from the DEG
results and queries the Enrichr API against configurable gene set libraries.
Supports upregulated, downregulated, or both directions independently — which
is critical for cancer data where tumour suppressors and activated oncogenes
represent biologically distinct pathway perturbations.

**Why it exists**

DEG results answer "which genes are different?" — GSEA answers "which biological
processes are different?" Pathway enrichment is the standard next step for
biological interpretation, connecting gene lists to the literature automatically.
Running up and down independently (rather than mixing both into one query)
produces interpretable results: an "up" enrichment means this cell type
*activates* this pathway, a "down" enrichment means it *suppresses* it.

**Steps performed (in order)**

1. Validate `direction` parameter and `gene_sets` names against Enrichr (warns, never crashes)
2. Build gene universe from `adata.var_names` (all detected genes = background for Fisher test)
3. For each group × direction: filter DEGs by `min_logfc`, `max_pval_adj`, `exclude_gene_prefixes` → build query gene list
4. Skip groups with < `min_genes` query genes — warns, stores in `skipped`
5. Run `gseapy.enrichr` per group × direction — derives Overlap count from Genes column (gseapy ≥1.0 dropped the Overlap column)
6. Tag each result row with `Direction` column
7. Build `summary_df` with top 3 terms per group (+ direction column when `direction="both"`)
8. Write provenance to `uns['omicsage_gsea']`

**`direction` parameter**

| Value | Query genes used | Result keys | Use case |
|-------|-----------------|-------------|----------|
| `"up"` (default) | logfc ≥ +min_logfc | `{group}` | Cell type marker enrichment |
| `"down"` | logfc ≤ -min_logfc | `{group}` | Suppressed pathways, tumour suppressor loss |
| `"both"` | Both, independently | `{group}__up`, `{group}__down` | Complete pathway picture for cancer/disease data |

**Key design decision — gene universe**

`exclude_gene_prefixes` filters the *query list only*, never the gene universe
(background). Removing genes from the background would artificially inflate
enrichment scores for pathways containing those genes. This is the statistically
correct approach for ORA.

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | AnnData from `deg()` — `var_names` used as gene universe |
| `deg_dict` | dict | required | `deg_dict` from `deg()` — must contain `'results'` key |
| `gene_sets` | list[str] | GO BP 2023, KEGG 2021, Reactome 2022 | Enrichr library names — any valid Enrichr name accepted |
| `min_logfc` | float | `0.25` | Minimum absolute log₂FC for query genes |
| `max_pval_adj` | float | `0.05` | Maximum adj. p-value for query genes |
| `top_n_genes` | int | None | Use top N DEGs per direction; None = use all passing filters |
| `min_genes` | int | `5` | Skip group×direction if fewer query genes |
| `organism` | str | `"human"` | Enrichr organism — `"human"` \| `"mouse"` \| etc. |
| `direction` | str | `"up"` | `"up"` \| `"down"` \| `"both"` — see table above |
| `exclude_gene_prefixes` | list[str] | None | Remove from query list only e.g. `["RPL", "RPS", "MT-"]` |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**Output**

```
adata_gsea.uns['omicsage_gsea']  →  provenance (includes direction, gene_sets, organism)

gsea_dict['results']    : dict   →  {group: DataFrame} for direction="up"/"down"
                                     {group__up: DataFrame, group__down: DataFrame} for direction="both"
                                     Each DataFrame: Term, Overlap, P-value, Adjusted P-value, Genes, Direction
gsea_dict['summary_df'] : DataFrame  →  top 3 terms per group; gains 'direction' column when direction="both"
gsea_dict['provenance'] : dict       →  same as uns['omicsage_gsea']
gsea_dict['skipped']    : list       →  (group, direction) tuples skipped due to < min_genes DEGs
```

**Usage**

```python
from pipeline.modules.downstream.gsea import gsea

# Standard immune data — upregulated pathways per cell type
adata_gsea, gsea_dict = gsea(
    adata_deg,
    deg_dict=deg_dict,
    gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"],
    direction="up",
    exclude_gene_prefixes=["RPL", "RPS", "MT-"],
    inplace=False,
)

# Cancer data — full picture (activated + suppressed pathways)
adata_gsea, gsea_dict = gsea(
    adata_deg,
    deg_dict=deg_dict,
    direction="both",
    gene_sets=["KEGG_2021_Human", "Reactome_2022", "DisGeNET"],
)

# Mouse dataset
adata_gsea, gsea_dict = gsea(
    adata_deg,
    deg_dict=deg_dict,
    gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Mouse"],
    organism="mouse",
)

# Inspect top pathways for T cells
print(gsea_dict["results"]["Tcm/Naive helper T cells"].head(5))

# For direction="both"
print(gsea_dict["results"]["Classical monocytes__up"].head(5))
print(gsea_dict["results"]["Classical monocytes__down"].head(5))
```

**Connects to**: `gsea_report.py` for the report

---

## 16. `gsea_report.py`

**What it does**

Generates a self-contained HTML report from `gsea()` output. All plots are
embedded as base64 PNGs. Handles both single-direction and `direction="both"`
results automatically — direction badges appear in tables and the summary when
running in "both" mode.

**Report contents**

| Section | What it shows |
|---------|--------------|
| Run summary | Stat cards (groups, gene sets, organism, sig. pathways) + query direction label + per-group sig. pathway counts table with direction badges (when applicable) |
| Top pathways table | Top 5 pathways per group; Genes Matched count, adj. p-value, gene list; direction badge in group cell when direction="both" |
| Bar charts | Top 10 pathways per group — horizontal bars sorted by −log₁₀(adj. p-value) |
| Bubble plot | Pathway × group matrix; size = genes matched, colour = −log₁₀(adj. p-value); when n_groups > max_bubble_groups, selects top N by sig. pathway count (not hard skip); excluded groups listed in visible note |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_dict` | dict | required | `gsea_dict` returned by `gsea()` |
| `output_path` | str | `"reports/gsea_report.html"` | HTML output path |
| `top_n_table` | int | 5 | Pathways per group in top pathways table |
| `top_n_bar` | int | 10 | Pathways per group in bar charts |
| `max_bubble_groups` | int | 9 | Max groups in bubble plot — excess handled by top-N selection |

**Usage**

```python
from reports.gsea_report import generate_gsea_report

report_path = generate_gsea_report(
    gsea_dict=gsea_dict,
    output_path="reports/gsea_report.html",
)
print(f"Report → {report_path}")
```

**Connects to**: (end of Phase 1 core pipeline) → batch correction or downstream analysis

---


## 17. `harmony_correct.py`

**What it does**

Runs Harmony batch integration on the PCA embedding produced by `reduce.py`.
Iteratively adjusts cell coordinates so that cells from different batches
intermix correctly while preserving biological variation. Recomputes the
neighbor graph and UMAP on the corrected embedding. Preserves the original
UMAP under a new key so before/after comparisons are always possible.

**Why it exists**

When a dataset spans multiple donors, sites, or sequencing runs, technical
batch effects can dominate the first principal components, causing cells to
cluster by batch rather than biology. Harmony is the standard, fast,
scalable correction method for scRNA-seq — it operates on the PCA embedding
(not the raw counts), is reversible, and integrates naturally into the
scanpy/AnnData workflow.

**Steps performed (in order)**

1. Validate inputs — checks `obsm['X_pca']` exists and `obs[batch_key]` has ≥ 2 unique values
2. Cap `n_pcs` at available PCA dimensions; logs a warning if capping occurs
3. Run `harmonypy.run_harmony` on the first `n_pcs` components of `X_pca`
4. Store corrected embedding in `obsm['X_pca_harmony']` — shape (n_cells × n_pcs)
5. Preserve existing `obsm['X_umap']` → `obsm['X_umap_precorrection']` (if present)
6. Recompute neighbor graph on corrected embedding → `uns['neighbors_harmony']`,
   `obsp['neighbors_harmony_connectivities']`, `obsp['neighbors_harmony_distances']`
7. Recompute UMAP on corrected graph → `obsm['X_umap_harmony']`
8. Write provenance to `uns['omicsage_harmony']`

**obsm / obsp layout after harmony_correct()**

| Slot | Contents |
|------|----------|
| `obsm['X_pca']` | Original PCA — **unchanged** |
| `obsm['X_pca_harmony']` | Harmony-corrected embedding (n_cells × n_pcs) |
| `obsm['X_umap_precorrection']` | Original UMAP preserved for comparison |
| `obsm['X_umap_harmony']` | UMAP recomputed on corrected embedding (n_cells × 2) |
| `obsp['neighbors_harmony_connectivities']` | kNN connectivity matrix on corrected embedding |
| `obsp['neighbors_harmony_distances']` | kNN distance matrix on corrected embedding |
| `uns['omicsage_harmony']` | Full provenance record |

**Important: `X_umap` is NOT overwritten.** The original UMAP from `reduce.py`
is preserved as `X_umap_precorrection`. The new UMAP is stored as `X_umap_harmony`.
This is a deliberate design decision — silently overwriting `X_umap` would make
before/after comparisons impossible and break any module that expects the
pre-correction UMAP.

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Reduced AnnData — must have `obsm['X_pca']` and `obs[batch_key]` |
| `batch_key` | str | `"batch"` | `obs` column encoding the batch variable |
| `n_pcs` | int | 50 | PCs to pass to Harmony — capped at available dimensions |
| `n_neighbors` | int | 15 | *k* for neighbor graph on corrected embedding |
| `umap_min_dist` | float | 0.3 | `min_dist` for UMAP |
| `random_state` | int | 0 | Reproducibility seed |
| `max_iter_harmony` | int | 50 | Maximum Harmony iterations |
| `theta` | float | 2.0 | Harmony diversity penalty — higher = stronger correction |
| `copy` | bool | False | If True, return a copy; otherwise modify in place |

**Usage**

```python
from pipeline.modules.integration.harmony_correct import harmony_correct

# Standard run on GSE194122 (10 donors = 10 batches)
adata = harmony_correct(
    adata,
    batch_key="batch",
    n_pcs=50,
    n_neighbors=15,
    umap_min_dist=0.3,
    random_state=0,
)

# Access outputs
print(adata.obsm["X_pca_harmony"].shape)      # (n_cells, 50)
print(adata.obsm["X_umap_harmony"].shape)     # (n_cells, 2)
print(adata.obsm["X_umap_precorrection"].shape) # (n_cells, 2) — original UMAP

# Provenance
import json
print(json.dumps(adata.uns["omicsage_harmony"], indent=2))
```

**Known implementation notes**

- `obs[batch_key]` is cast to `str` internally — safe when the column is a
  pandas `Categorical` (avoids `NotImplementedError: isna is not defined for MultiIndex`)
- `ho.Z_corr` from `harmonypy` is already `(n_cells, n_pcs)` — do **not** transpose;
  earlier versions of this module incorrectly used `.T` which caused shape errors
- `neighbors_key="neighbors_harmony"` keeps the uncorrected graph intact in
  `obsp['connectivities']` — allows direct before/after comparison
- Dependency: `pip install harmonypy`

**Connects to**: `harmony_report.py` for the report, then clustering on `neighbors_harmony`

---

## 18. `harmony_report.py`

**What it does**

Generates a self-contained HTML report from `harmony_correct()` output.
Verifies output keys, visualises batch mixing before and after correction,
computes a quantitative mixing score, and profiles how much each PC was
shifted by the correction. All plots are embedded as base64 PNGs.
Matches the style of `gsea_report.py`.

**Report contents**

| Section | What it shows |
|---------|--------------|
| Run summary | Stat cards (cells, genes, batches, PCs corrected, k, elapsed) + output key verification table (✓ / ✗ MISSING) + run parameters (`batch_key`, `theta`, `max_iter`, `umap_min_dist`) |
| Batch composition | Horizontal bar chart of cells per batch + table with cell count and % of total |
| UMAP embeddings | Side-by-side: raw PCA PC1–PC2 (before, coloured by batch) vs corrected UMAP (after, coloured by batch); second panel: UMAP coloured by Harmony PC1 value (shows correction depth) |
| Batch mixing metrics | Histogram of per-cell same-batch neighbour fraction; mean / median / expected (= 1/n_batches) stats; normalised mixing score (0–1) with colour-coded interpretation note |
| Per-PC correction shift | Bar chart of mean `|X_pca − X_pca_harmony|` per PC; top 5 most-shifted PCs table |

**Mixing score interpretation**

| Score | Meaning |
|-------|---------|
| ≥ 0.8 | ✓ Batches are well integrated |
| 0.5–0.8 | ⚠ Integration is moderate — consider increasing `theta` |
| < 0.5 | ✗ Batches are poorly integrated — check `batch_key` or try higher `theta` |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | AnnData returned by `harmony_correct()` |
| `output_path` | str | `"reports/harmony_report.html"` | HTML output path |
| `top_n_pcs` | int | 20 | PCs to show in the per-PC shift chart |
| `max_batches_bar` | int | 20 | Max batches labelled individually in the composition bar chart |

**Usage**

```python
from reports.harmony_report import generate_harmony_report

report_path = generate_harmony_report(
    adata=adata,
    output_path="reports/harmony_report.html",
)
print(f"Report → {report_path}")
```

**Connects to**: clustering on `neighbors_harmony` graph — next step in pipeline

---

## 19. `pseudobulk_deg.py`

**What it does**

Aggregates raw counts per (cell_type, donor) into bulk-like pseudo-samples and
runs DESeq2 Wald tests via `pydeseq2`, one-vs-rest per cell type. Complements
the Wilcoxon results from `deg.py` — DESeq2 is more conservative but correctly
controls for donor-level variation, making it the preferred method for publication
figures and datasets with replicated donors.

**Why it exists**

The Wilcoxon test in `deg.py` treats every cell as an independent observation.
With thousands of cells per cell type, this inflates statistical power dramatically
and produces many false positives driven by cell number rather than biology
(the "pseudo-replication" problem). Pseudobulk DEG aggregates cells from the same
donor into a single sample, restoring biological replication as the unit of
inference. DESeq2's negative binomial model then handles count noise correctly.

**Steps performed (in order)**

1. Validate obs columns (`groupby`, `donor_key`) and layers key (`counts_layer`)
2. Cast `obs[donor_key]` and `obs[groupby]` to `str` — avoids Categorical/MultiIndex errors
3. Extract raw counts matrix from `layers[counts_layer]` as dense float64
4. Aggregate: sum counts per `(cell_type, donor)` combo; drop combos below `min_cells`
5. For each cell type: check `n_target` ≥ `min_samples` and `n_rest` ≥ `min_samples`; skip with `UserWarning` if not
6. Build `DeseqDataSet` with `condition` covariate (`target_ct` vs `"rest"`)
7. Run `dds.deseq2()` + `DeseqStats(contrast=["condition", target_ct, "rest"])`
8. Rename pydeseq2 columns → deg.py schema: `log2FoldChange→logfc`, `pvalue→pval`, `padj→pval_adj`, `stat→score`
9. Drop NaN `padj` rows (low-count outliers excluded by DESeq2)
10. Apply `min_logfc` / `max_pval_adj` thresholds; apply `exclude_gene_prefixes`
11. Build `summary_df` (top 5 per group, long format) — identical schema to `deg.py`
12. Write provenance to `uns['omicsage_pseudobulk_deg']`

**Output schema — identical to `deg.py` deg_dict**

| Key | Type | Contents |
|-----|------|----------|
| `results` | dict | `{group: DataFrame(gene, score, pval, logfc, pval_adj)}` — significant DEGs only |
| `summary_df` | DataFrame | Long-format top 5 DEGs per group: group, rank, gene, logfc, pval_adj |
| `provenance` | dict | Same metadata as `uns['omicsage_pseudobulk_deg']` |
| `pairwise` | dict | `{}` — not implemented; kept for `deg_report.py` compatibility |
| `skipped` | dict | `{group: reason}` — cell types skipped due to too few donor pseudo-samples |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Must have `obs[groupby]`, `obs[donor_key]`, `layers[counts_layer]` |
| `groupby` | str | `"cell_type_vote"` | obs column for cell-type grouping |
| `donor_key` | str | `"batch"` | obs column for donor identity |
| `counts_layer` | str | `"counts"` | Layer with raw integer counts — do **not** pass `"logcounts"` |
| `min_cells` | int | 10 | Minimum cells per (cell_type, donor) pseudo-sample |
| `min_samples` | int | 3 | Minimum donor pseudo-samples per cell type (both target and rest groups) |
| `min_logfc` | float | `0.25` | Minimum absolute log₂FC threshold |
| `max_pval_adj` | float | `0.05` | Maximum BH-corrected adjusted p-value (DESeq2 padj) |
| `n_top_genes` | int | 500 | Maximum significant genes to retain per group after thresholding |
| `exclude_gene_prefixes` | list[str] | None | Gene prefixes to remove post-threshold e.g. `["RPL", "RPS", "MT-"]` |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**Usage**

```python
from pipeline.modules.downstream.pseudobulk_deg import pseudobulk_deg

adata_pb, pb_dict = pseudobulk_deg(
    adata_annotated,
    groupby="cell_type_vote",
    donor_key="batch",
    counts_layer="counts",
    min_cells=10,
    min_samples=3,
    min_logfc=0.25,
    max_pval_adj=0.05,
    inplace=False,
)

# Check which groups were skipped
print(pb_dict["skipped"])

# Inspect significant DEGs for T cells
print(pb_dict["results"]["T_cell"].head(10))

# Top 5 DEGs per group
print(pb_dict["summary_df"])
```

**Critical implementation notes**

- `counts_layer` must hold raw **integer** counts — DESeq2 requires integers. Passing `"logcounts"` will produce silently wrong results because summing log-normalised values is not a valid bulk aggregation.
- `obs[donor_key]` is cast to `str` internally — same Categorical/MultiIndex bug that affects Harmony and GSEA.
- When `target_ct` is the target group, the rest group includes *all other cell types*. Both `n_target` and `n_rest` are checked independently against `min_samples`. If either fails, the cell type is skipped (not crashed).
- pydeseq2 drops genes with zero total counts per comparison — the gene list in each result DataFrame may therefore differ slightly from the full `var_names`.
- `uns['rank_genes_groups']` is NOT populated — pseudobulk DEG does not use scanpy's `rank_genes_groups`. The dot plot in `pseudobulk_deg_report.py` uses `layers['logcounts']` for display instead.
- Dependency: `pip install pydeseq2`

**Connects to**: `pseudobulk_deg_report.py` for the report

---

## 20. `pseudobulk_deg_report.py`

**What it does**

Generates a self-contained HTML report from `pseudobulk_deg()` output. Identical
CSS and page skeleton to `deg_report.py` — same dark-blue header, same stat cards,
same volcano grid — with two additional sections not present in `deg_report.py`:
a **Skipped Groups** table explaining why each cell type was excluded, and
pseudobulk-specific stat cards (`donor_key`, `counts_layer`, `min_cells`, `min_samples`).

Does **not** require or use `uns['rank_genes_groups']` — reads
`uns['omicsage_pseudobulk_deg']` natively. The dot plot is rendered from
`layers['logcounts']` directly.

**Report contents**

| Section | What it shows |
|---------|--------------|
| Run summary | Stat cards: groups tested, groups skipped, method, groupby, donor key, counts layer, min cells/samples, thresholds, total significant DEGs; per-group significant DEG counts table; prefix exclusion note (if applicable); pseudobulk method explanation note |
| Skipped groups | Table of cell types that failed `min_samples` check + exact reason; "none skipped" confirmation when all passed |
| Top DEGs per group | Rowspan table — top 5 significant genes per group, ranked by DESeq2 padj; Direction column (▲ Up / ▼ Down) + log₂FC coloured red/blue |
| Volcano plots | One per group (default max 20, sorted by DEG count); same style as `deg_report.py` |
| Dot plot | `sc.pl.dotplot` using `layers['logcounts']` for display — top 5 DEGs per group |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | AnnData returned by `pseudobulk_deg()` — must have `uns['omicsage_pseudobulk_deg']` |
| `pb_dict` | dict | required | `pb_dict` returned by `pseudobulk_deg()` |
| `output_path` | str | `"reports/pseudobulk_deg_report.html"` | HTML output path |
| `top_n_volcano` | int | 10 | Genes to label on each volcano plot |
| `top_n_dotplot` | int | 5 | Top DEGs per group in dot plot |
| `max_volcano_groups` | int | 20 | Max volcano plots; excess sorted by DEG count |

**Usage**

```python
from reports.pseudobulk_deg_report import generate_pseudobulk_deg_report

report_path = generate_pseudobulk_deg_report(
    adata=adata_pb,
    pb_dict=pb_dict,
    output_path="reports/pseudobulk_deg_report.html",
)
print(f"Report → {report_path}")
```

**Connects to**: (end of Phase 1 DEG branch) → HCC milestone benchmark

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
        │       │       normalization_report.py → reports/normalization_report.html
        │       │
        │       ▼
        │   reduce.py     → obsm['X_pca'] + obsm['X_umap'] + obsp['connectivities']
        │       │                   │
        │       │                   ▼
        │       │       reduce_report.py → reports/reduce_report.html
        │       │
        │       ▼
        │   cluster.py    → obs['leiden_*'] + obs['leiden']
        │       │                   │
        │       │                   ▼
        │       │       cluster_report.py → reports/cluster_report.html
        │       │
        │       ▼
        │   annotate.py   → obs['cell_type_vote'] + obs['cell_type_confidence']
        │       │                   │
        │       │                   ▼
        │       │       annotate_report.py → reports/annotate_report.html
        │       │
        │       ▼
        │   deg.py        → deg_dict['results'] + uns['rank_genes_groups']
        │       │           (rankby_abs=True — both directions)
        │       │                   │
        │       │                   ▼
        │       │       deg_report.py → reports/deg_report.html
        │       │
        │       ▼
        │   gsea.py       → gsea_dict['results'] (direction: up | down | both)
        │       │                   │
        │       │                   ▼
        │       │       gsea_report.py → reports/gsea_report.html
        │       │
        │       ▼
        │   harmony_correct.py  → obsm['X_pca_harmony'] + obsm['X_umap_harmony']
        │       │                     obsm['X_umap_precorrection'] preserved
        │       │                   │
        │       │                   ▼
        │       │       harmony_report.py → reports/harmony_report.html
        │       │
        │       ▼
        │   cluster.py (neighbors_key='neighbors_harmony')
        │       │         → obs['leiden_harmony'] + obs['leiden_harmony_<res>']
        │       │         → compute_ari('leiden', 'leiden_harmony')
        │       │
        │       ▼
        │   pseudobulk_deg.py   → pb_dict['results'] (DESeq2 Wald, one-vs-rest)
        │       │                   layers['counts'] aggregated per (cell_type, donor)
        │       │                   uns['omicsage_pseudobulk_deg']
        │       │                   │
        │       │                   ▼
        │       │       pseudobulk_deg_report.py → reports/pseudobulk_deg_report.html
        │       │
        │       ▼
        │   [MILESTONE: Wang et al. 2025 HCC benchmark]  ← NEXT
        │
        ├── mdata["adt"]  → ADT QC + CLR normalization (future phase)
        └── mdata["atac"] → ATAC QC (Phase 4)
```

---

## Tests

| Test file | Module | Tests | What it covers |
|-----------|--------|-------|---------------|
| `tests/test_phase0_structure.py` | — | — | Repo structure, imports, config schema |
| `tests/test_ingest.py` | `ingest.py` | — | Format detection, raw count extraction, all three loaders |
| `tests/test_qc.py` | `qc.py` | 42 | MT detection, metrics, filtering, Scrublet, ground-truth validation, modality detection, MuData structure |
| `tests/test_normalize.py` | `normalize.py` | 12 | Raw count preservation, normalization correctness, log1p, HVG selection, batch_key, logcounts layer, provenance, mutation guard |
| `tests/test_reduce.py` | `reduce.py` | 12 | PCA/UMAP shapes, HVG-only PCA, neighbor graph, provenance, inplace guard, t-SNE optional, PC selection methods |
| `tests/test_cluster.py` | `cluster.py` | 16 | Leiden labels, all resolutions computed, n_clusters bounds, silhouette scores, best resolution selection, provenance (incl. neighbors_key/cluster_key), inplace guard, neighbors_key harmony routing, cluster_key isolation, compute_ari, missing harmony connectivities error |
| `tests/test_annotate.py` | `annotate.py` | 18 (+1 skip) | obs columns written, provenance keys, confidence range, marker scoring, vote consensus, inplace guard, ground-truth preservation |
| `tests/test_deg.py` | `deg.py` | 11 | Return types, provenance keys, column names, pval range, per-group results, threshold filtering, inplace guard, small-group warning, leiden fallback |
| `tests/test_gsea.py` | `gsea.py` | 8 | Return types, provenance keys, output columns, pval range, group accounting, min_genes skip + warning, inplace guard, string gene_sets param — all Enrichr calls mocked (CI-safe) |
| `tests/test_harmony.py` | `harmony_correct.py` | 13 | Embedding created, shape (n_cells × n_pcs), UMAP recomputed as X_umap_harmony, X_umap_precorrection preserved and unchanged, neighbors key stored, provenance keys + values, in-place modification, n_pcs capping, missing X_pca error, missing batch_key error, single-batch error, custom batch_key, two-batch correction |
| `tests/test_pseudobulk_deg.py` | `pseudobulk_deg.py` | 14 | Return types, provenance keys, pb_dict schema matches deg_dict (results/summary_df/provenance/pairwise/skipped), result DataFrame columns, pval range, all groups accounted for (results ∪ skipped), threshold filtering (logfc + pval_adj), inplace guard, cell type with too few donors skipped + UserWarning issued, missing counts_layer error, missing groupby column error |

Run all tests:

```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~231 passed, 1–2 skipped
```

Run a single module's tests:

```bash
python -m pytest tests/test_cluster.py -v           # 16 passed
python -m pytest tests/test_harmony.py -v           # 13 passed
python -m pytest tests/test_deg.py -v               # 11 passed
python -m pytest tests/test_annotate.py -v          # 18 passed, 1 skipped
python -m pytest tests/test_pseudobulk_deg.py -v    # 14 passed
```

---

## 21. `run_scrna_pipeline.py` + `config/runs/`

**What it does**

Generic pipeline runner. Executes any subset of the 10 Phase 1 steps for any
dataset. All dataset identity, paths, and parameters live in a per-dataset YAML
config file — the runner itself never needs to be edited.

**Location**: `run_scrna_pipeline.py` (repo root)

**Config files**: `config/runs/<dataset_id>.yaml` — one file per dataset

---

### Config file structure

```yaml
dataset:
  id: GSE194122
  name: "BMMC CITE-seq (NeurIPS 2021)"
  modality: cite        # cite | scrna | atac | multiome
  organism: human

paths:
  raw_input: data/benchmark/GSE194122_cite_raw_only.h5ad
  processed_dir: data/processed/GSE194122
  reports_dir: reports/GSE194122

steps:
  qc:
    enabled: true
    params:
      min_genes: 200
      max_genes: 2500
      max_mt_pct: 5.0
      remove_doublets: true

  normalize:
    enabled: true
    params:
      batch_key: batch
      target_sum: 10000
      n_top_genes: 2000
      hvg_flavor: seurat

  # ... one block per step, same pattern
  harmony:
    enabled: false        # set false when no batch effects
```

Every step block supports an `input_override` field to inject a pre-existing
file as the step's input, bypassing predecessor output:

```yaml
reduce:
  enabled: true
  input_override: data/external/already_normalized.h5ad
```

Full parameter reference: `config/schema.yaml`

---

### Step output files

| Step | File saved | Report saved |
|------|-----------|-------------|
| qc | `01_qc.h5ad` + `01_qc_adt.h5ad` | `01_qc_report.html` |
| normalize | `02_normalized.h5ad` | `02_normalization_report.html` |
| reduce | `03_reduced.h5ad` | `03_reduce_report.html` |
| cluster | `04_clustered.h5ad` | `04_cluster_report.html` |
| annotate | `05_annotated.h5ad` | `05_annotate_report.html` |
| deg | `06_deg.h5ad` | `06_deg_report.html` |
| gsea | `07_gsea.h5ad` | `07_gsea_report.html` |
| harmony | `08_harmony.h5ad` | `08_harmony_report.html` |
| cluster_harmony | `09_harmony_clustered.h5ad` | — |
| pseudobulk | `10_pseudobulk_deg.h5ad` | `10_pseudobulk_deg_report.html` |

All paths are under `paths.processed_dir` and `paths.reports_dir` from the config.

`pseudobulk` reads from `annotate` output (step 5), not `cluster_harmony` (step 9) —
it needs `obs['cell_type_vote']`, `obs['batch']`, and `layers['counts']` which live
on the annotated file. Harmony does not add or modify these columns.

---

### CLI reference

```bash
# Full pipeline
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml

# Stop at a checkpoint (inclusive) — then inspect the report
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --to-step cluster

# Resume from a checkpoint (inclusive)
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --from-step annotate

# Run a specific range
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --from-step normalize --to-step reduce

# Run exactly one step
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --step normalize

# Valid step names (in order):
# qc  normalize  reduce  cluster  annotate  deg  gsea  harmony  cluster_harmony  pseudobulk
```

---

### Interactive checkpoint workflow (cluster resolution)

The standard single-cell workflow requires inspecting clusters before continuing.
The runner supports this explicitly:

```bash
# 1. Run through clustering
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --to-step cluster

# 2. Open reports/GSE194122/04_cluster_report.html
#    Decide you want resolution 0.8 instead of the auto-selected value

# 3. Add to config under cluster.params:
#      resolution_override: 0.8

# 4. Re-run clustering with the chosen resolution
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --step cluster

# 5. Continue from annotation onward
python run_scrna_pipeline.py --config config/runs/GSE194122.yaml --from-step annotate
```

---

### Common config patterns

**No batch effects — skip Harmony:**
```yaml
harmony:         { enabled: false }
cluster_harmony: { enabled: false }
```

**Pre-processed data — start from reduce:**
```yaml
qc:        { enabled: false }
normalize: { enabled: false }
reduce:
  enabled: true
  input_override: data/external/already_normalized.h5ad
```

**New dataset:** copy `config/runs/GSE194122.yaml` → `config/runs/GSE166635.yaml`,
update `dataset.id`, `paths.raw_input`, and any parameters that differ.
Zero Python required.

---

### Validation

The runner validates all inputs before executing any step. If a step's required
input is missing and no `input_override` is provided, it exits with a clear
message before touching any data:

```
[OmicSage] Validation failed — fix these before running:

  ✗  [reduce] requires output of 'normalize' at data/processed/GSE194122/02_normalized.h5ad
       Options:
         • Run 'normalize' first
         • Add 'input_override' under steps.reduce in your config
```

---

### Caching

If a step's output `.h5ad` already exists it is skipped (cached) and the
existing path is passed to the next step. The one exception is `cluster` when
`resolution_override` is set — in that case the step always re-runs to apply
the analyst's chosen resolution.

---

### Logging

```bat
@echo off
set PYTHONUTF8=1
echo. >> logs\GSE194122.log
echo ======================================================== >> logs\GSE194122.log
echo RUN %date% %time% >> logs\GSE194122.log
echo ======================================================== >> logs\GSE194122.log
python run_scrna_pipeline.py --config configs\GSE194122.yaml >> logs\GSE194122.log 2>&1
```

The runner prints start time, end time, and elapsed time (`HHh MMm SSs`) at
the end of every run.

---

## 22. `reports/combined_report.py`

**What it does**

Assembles all individual step HTML reports into a single self-contained
tabbed HTML file. Each pipeline step gets one tab. Only tabs for reports
that actually exist on disk are shown — partial runs produce a partial
combined report without errors.

**Why it exists**

After a full pipeline run, the `reports/<dataset>/` directory contains
up to 9 separate HTML files. A biologist or collaborator has no obvious
entry point and must open each file individually. `combined_report.py`
provides a single file that is the complete analysis record, with a
progress bar showing how much of the pipeline has been run.

**Design decisions**

- Zero changes to existing report generators — the combiner reads their
  output after the fact rather than modifying how they write
- Extracts `<main>` content from each step report; falls back to `<body>`
  minus header/footer if no `<main>` tag is present
- Step-specific CSS (volcano grid sizing, dot plot widths, etc.) is
  preserved by extracting and re-embedding `<style>` blocks
- No new dependencies — uses only Python stdlib (`re`, `datetime`, `pathlib`)

**Tab registry** (in display order)

| Filename | Tab label |
|----------|-----------|
| `01_qc_report.html` | QC |
| `02_normalization_report.html` | Normalize |
| `03_reduce_report.html` | Reduce |
| `04_cluster_report.html` | Cluster |
| `05_annotate_report.html` | Annotate |
| `06_deg_report.html` | DEG |
| `07_gsea_report.html` | GSEA |
| `08_harmony_report.html` | Harmony |
| `10_pseudobulk_deg_report.html` | Pseudobulk |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reports_dir` | Path | required | Directory containing step HTML reports |
| `dataset_name` | str | `"OmicSage Analysis"` | Shown in report header |
| `output_path` | Path | `reports_dir / "00_combined_report.html"` | Output path |

**Output**

```
00_combined_report.html   →  single self-contained tabbed HTML file
                              sorts to top of folder due to "00_" prefix
```

**Header features**

- Dataset name + generation timestamp
- Progress bar: `N of 9 pipeline steps complete` + percentage
- Step labels listed inline

**Navigation features**

- Click any tab to switch panels
- Left/right arrow keys navigate between tabs
- Tab button scrolls into view on mobile (many tabs case)

**Usage — called automatically by run_scrna_pipeline.py**

```python
from reports.combined_report import generate_combined_report

generate_combined_report(
    reports_dir=Path("reports/GSE166635"),
    dataset_name="GSE166635 — HCC",
)
# → reports/GSE166635/00_combined_report.html
```

**Usage — standalone rebuild from existing reports**

```bash
python -m reports.combined_report \
  --reports-dir reports/GSE166635 \
  --dataset-name "GSE166635 — HCC"

# Custom output path:
python -m reports.combined_report \
  --reports-dir reports/GSE166635 \
  --dataset-name "GSE166635 — HCC" \
  --output reports/GSE166635/combined.html
```

**Wiring into run_scrna_pipeline.py**

Add at the end of `main()`, just before the footer print block:

```python
from reports.combined_report import generate_combined_report
generate_combined_report(
    reports_dir=reports_dir,
    dataset_name=dataset_name,
    output_path=reports_dir / "00_combined_report.html",
)
```

**Connects to**: all step report generators — reads their output files,
produces no `.h5ad` output of its own.
