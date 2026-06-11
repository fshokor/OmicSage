# OmicSage — Spatial Transcriptomics Module Documentation
> Phase 7 — Spatial Transcriptomics (Sessions 1–5) + Session B extensions (Session B)
> Last updated: June 2026

This document covers every script produced in Phase 7 and the Session B
extension: pipeline modules, report generators, the combined report, the runner,
the config, the tests, and the benchmark download script.

---

## Table of Contents

1. [Pipeline Architecture](#1-pipeline-architecture)
2. [Pipeline Modules](#2-pipeline-modules)
   - [spatial_ingest.py](#21-spatial_ingestpy)
   - [spatial_qc.py](#22-spatial_qcpy)
   - [spatial_reduce.py](#23-spatial_reducepy)
   - [spatial_cluster.py](#24-spatial_clusterpy)
   - [spatial_deconvolve.py](#25-spatial_deconvolvepy)
   - [spatial_downstream.py](#26-spatial_downstreampy)
   - [spatial_impute.py](#27-spatial_imputepy)
3. [Report Templates](#3-report-templates)
   - [spatial_qc_report.py](#31-spatial_qc_reportpy)
   - [spatial_reduce_report.py](#32-spatial_reduce_reportpy)
   - [spatial_cluster_report.py](#33-spatial_cluster_reportpy)
   - [spatial_deconvolve_report.py](#34-spatial_deconvolve_reportpy)
   - [spatial_downstream_report.py](#35-spatial_downstream_reportpy)
   - [spatial_impute_report.py](#36-spatial_impute_reportpy)
4. [Combined Report](#4-combined-report)
5. [Pipeline Runner](#5-pipeline-runner)
6. [Config Files](#6-config-files)
7. [Tests](#7-tests)
8. [Benchmark Data Download](#8-benchmark-data-download)
9. [AnnData State Reference](#9-anndata-state-reference)
10. [Key Decisions](#10-key-decisions)

---

## 1. Pipeline Architecture

```
spatial_ingest.py
      │  standard AnnData: obsm["spatial"], uns["spatial"], layers["counts"]
      │  Supported: visium | h5ad | benchmark | visium_hd | xenium
      ▼
spatial_qc.py
      │  obs QC columns, qc_pass mask, filtered spots
      ▼
spatial_reduce.py
      │  X normalized, var["highly_variable"], obsm["X_pca"],
      │  obsp["spatial_connectivities"]
      ▼
spatial_cluster.py
      │  obs["spatial_cluster"], uns["moranI"]
      ▼
spatial_deconvolve.py
      │  obsm["q05_cell_abundance_w_sf"], obs[cell_type] columns
      │  (graceful skip when no reference — downstream still runs)
      ▼
spatial_downstream.py
      │  obs["region_cluster"], uns["celltype_marker_genes"],
      │  uns["celltype_svg"], uns["*_nhood_enrichment"],
      │  uns["*_co_occurrence"], uns["*_ligrec"], uns["svg_gsea"]
      ▼
spatial_impute.py          ← Session B addition
      │  obsm["imputed_expression"] (float32 array, spots × n_genes)
      │  uns["omicsage_spatial_impute"] with genes_imputed list
      │  (graceful skip when sc_reference_path absent or enabled=false)
```

Every module follows the same contract:
- Returns `(AnnData, dict)` — the modified AnnData and a provenance dict
- `inplace=False` by default — operates on a copy
- Writes provenance to `adata.uns["omicsage_<module>"]`
- All parameters have sensible defaults drawn from the Kuppe benchmark

---

## 2. Pipeline Modules

### 2.1 `spatial_ingest.py`

**Location:** `pipeline/modules/spatial/spatial_ingest.py`
**Phase 7 Session:** 1 (updated Session 4, updated Session B)

**Purpose:** Unified entry point for loading spatial transcriptomics data
from any supported technology into a standard AnnData object. Technology
is specified via `spatial_type` or auto-detected from the source path.

**Public API:**

```python
spatial_ingest(
    source: str,                         # path or "benchmark"
    spatial_type: str = "auto",          # "auto" | "visium" | "h5ad" | "benchmark"
                                         # | "visium_hd" | "xenium" | "merfish" | "codex"
    counts_file: str = "filtered_feature_bc_matrix.h5",
    library_id: Optional[str] = None,
    library_key: Optional[str] = None,
    load_images: bool = True,
    bin_size: int = 8,                   # Visium HD only: 2 | 8 | 16 µm
    inplace: bool = False,
) -> tuple[AnnData, dict]

list_supported_types() -> dict[str, str]  # maps type → "implemented"|"planned"
```

**Supported technologies:**

| Type | Status | Source format | Loader |
|------|--------|---------------|--------|
| `visium` | ✅ Implemented | Space Ranger output directory | `sq.read.visium` |
| `h5ad` | ✅ Implemented | Pre-built `.h5ad` file | `sc.read_h5ad` |
| `benchmark` | ✅ Implemented | `squidpy.datasets.visium_hne_adata()` | squidpy datasets |
| `visium_hd` | ✅ Implemented | `binned_outputs/` directory | `spatialdata_io.visium_hd` |
| `xenium` | ✅ Implemented | Xenium output directory | `spatialdata_io.xenium` |
| `merfish` | 🔲 Planned | `cell_by_gene.csv` | — |
| `codex` | 🔲 Planned | `.csv` protein markers | — |

**Auto-detection fingerprints** (checked in order):
1. `source == "benchmark"` → benchmark
2. `source.endswith(".h5ad")` → h5ad
3. Directory contains `transcripts.parquet` → **xenium**
4. Directory contains `cell_by_gene.csv` → merfish
5. Directory contains `binned_outputs/` → **visium_hd**
6. Directory contains `spatial/` subfolder → visium
7. File ends with `.csv` → codex

**Visium HD loader (`_load_visium_hd`)** — *Session B*:
- Requires `pip install spatialdata-io`; raises `ImportError` with install
  instructions when absent
- Calls `spatialdata_io.visium_hd(path, bin_size=N, load_segmentations_only=False)`
- Table key is `f"square_{bin_size:03d}um"` (e.g. `"square_008um"`)
- `obsm["spatial"]` is already set by spatialdata-io from `pxl_col_in_fullres`
  / `pxl_row_in_fullres` — no extra extraction needed
- Builds minimal `uns["spatial"][lib_id]` stub so downstream squidpy tools pass
  validation; includes `scalefactors["spot_diameter_fullres"] = bin_size`
- `KeyError` raised with useful message when the requested bin_size table is absent

**Xenium loader (`_load_xenium`)** — *Session B*:
- Requires `pip install spatialdata-io`
- Calls `spatialdata_io.xenium(path, cells_table=True, cells_boundaries=False,
  nucleus_boundaries=False, transcripts=False, morphology_mip=False, ...)`
  — heavy assets (segmentation masks, morphology images, transcripts) are
  deliberately skipped by default to keep memory usage low
- Table key is always `"table"`
- `obsm["spatial"]` already set from `x_centroid` / `y_centroid`
- Builds minimal `uns["spatial"][lib_id]` stub with `"platform": "xenium"` in
  `metadata`
- Key difference from Visium: cell-level (not spot-level), targeted gene panel

**h5ad loading contract** (updated Session 4 for Kuppe data):
- If `layers["counts"]` absent → saves `X.copy()` to `layers["counts"]`
- If `var["gene_ids"]` present → swaps `var_names` to ENSEMBL IDs, preserves symbols in `var["feature_name"]`
- If any `MT-` genes present → strips them into `obsm["MT"]` and removes from `X`

**Output AnnData contract** (all implemented types):
```
obsm["spatial"]                             coordinates (n_obs, 2)
uns["spatial"][library_id]["images"]        tissue images (empty dict for Visium HD / Xenium by default)
uns["spatial"][library_id]["scalefactors"]  scale factors
layers["counts"]                            raw integer counts
uns["omicsage_spatial_ingest"]              provenance dict
```

**Provenance keys:** `source`, `spatial_type`, `technology_notes`, `counts_file`,
`library_id`, `library_key`, `load_images`, `bin_size`, `n_obs`, `n_vars`, `timestamp`

---

### 2.2 `spatial_qc.py`

**Location:** `pipeline/modules/spatial/spatial_qc.py`
**Phase 7 Session:** 1

**Purpose:** Compute QC metrics for Visium spots / Xenium cells and optionally
filter low-quality observations based on count/gene/MT% thresholds.

**Public API:**

```python
spatial_qc(
    adata: AnnData,
    min_counts: int = 500,
    max_counts: int = 100_000,
    min_genes: int = 200,
    max_genes: int = 10_000,
    max_mt_pct: float = 20.0,
    mt_prefix: str = "MT-",              # "MT-" human, "mt-" mouse
    filter_spots: bool = True,
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**Technology-specific threshold guidance:**
- **Visium:** `min_counts 500`, `min_genes 200` — standard
- **Visium HD 8µm:** `min_counts 50–200`, `min_genes 50` — smaller bins
- **Xenium:** `min_counts 5–50`, `min_genes 5–20` — panel-limited

**What it does:**
1. Annotates `var["mt"]` using `mt_prefix`
2. Runs `sc.pp.calculate_qc_metrics` → adds `total_counts`, `n_genes_by_counts`, `pct_counts_mt` to `obs`
3. Computes per-threshold failure counts (before filtering)
4. Builds `obs["qc_pass"]` boolean mask
5. Filters spots if `filter_spots=True`
6. Computes summary stats (mean/median/std/min/max) on retained spots

**Input requirements:** `obsm["spatial"]` must be present.

**Output obs columns added:**
```
total_counts          total UMI counts per spot/cell
n_genes_by_counts     number of detected genes per spot/cell
pct_counts_mt         mitochondrial gene percentage
qc_pass               bool — True if spot passes all thresholds
```

**Provenance keys:** all threshold parameters, `n_spots_before`, `n_spots_after`,
`n_spots_removed`, per-filter removal counts, summary stats per metric.

**Note on Kuppe Visium data:** After `spatial_ingest` strips MT genes into
`obsm["MT"]`, `pct_counts_mt` will be 0 for all spots. The `max_mt_pct`
threshold is vacuously satisfied — this is correct behaviour.

---

### 2.3 `spatial_reduce.py`

**Location:** `pipeline/modules/spatial/spatial_reduce.py`
**Phase 7 Session:** 2

**Purpose:** Normalize counts, select highly variable genes, compute PCA,
and build the spatial neighbours graph.

**Public API:**

```python
spatial_reduce(
    adata: AnnData,
    n_top_genes: int = 3000,
    n_comps: int = 50,
    n_neighbors: int = 6,                # 6 = Visium hex grid; use ~10 for Xenium
    coord_type: Optional[str] = None,    # None = auto-detect grid
    normalize_total: bool = True,
    target_sum: float = 1e4,
    log1p: bool = True,
    flavor: str = "seurat",
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**Technology-specific guidance:**
- **Visium / Visium HD:** `n_top_genes 3000`, `n_neighbors 6` (hex grid)
- **Xenium:** `n_top_genes` should equal or approach panel size; `n_neighbors 10`

**What it does:**
1. **Raw count restoration** — if `layers["counts"]` present, resets `X` to raw before normalizing
2. **Benchmark detection** — skips normalization when data already processed
3. **Normalization** — `sc.pp.normalize_total` → `sc.pp.log1p`
4. **HVG selection** — `sc.pp.highly_variable_genes(flavor="seurat", n_top_genes=...)`
5. **PCA** — `sc.pp.pca(use_highly_variable=True, svd_solver="arpack")`
6. **Spatial neighbours** — `sq.gr.spatial_neighbors(coord_type=None)` — auto-detects grid

**Output keys added:**
```
X                              normalized + log1p counts
layers["counts"]               raw counts (preserved)
var["highly_variable"]         HVG boolean mask
obsm["X_pca"]                  PCA embedding (n_obs × n_comps)
uns["pca"]["variance_ratio"]   explained variance per PC
obsp["spatial_connectivities"] spatial adjacency matrix
obsp["spatial_distances"]      spatial distances
uns["spatial_neighbors"]       spatial graph metadata
uns["omicsage_spatial_reduce"] provenance dict
```

**Provenance keys:** all parameters, `skipped_normalization`, `n_hvgs`,
`n_comps_computed`, `pca_variance_ratio_top10`, `pca_cumulative_variance_top10`,
`spatial_graph_n_edges`, `spatial_graph_mean_neighbors`.

---

### 2.4 `spatial_cluster.py`

**Location:** `pipeline/modules/spatial/spatial_cluster.py`
**Phase 7 Session:** 3

**Purpose:** Leiden clustering on transcriptomic similarity (KNN graph from
PCA), plus Moran's I spatially variable gene detection on the spatial graph.

**Public API:**

```python
spatial_cluster(
    adata: AnnData,
    resolution: float = 0.5,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    random_state: int = 0,
    cluster_key: str = "spatial_cluster",
    annotation_map: Optional[dict] = None,
    run_svg: bool = True,
    svg_n_genes: Optional[int] = None,
    svg_n_jobs: int = 1,
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**Key design decision:** Leiden uses the transcriptomic KNN graph, NOT the
spatial graph. The spatial graph is reserved for Moran's I.

**Output keys added:**
```
obs[cluster_key]                Leiden cluster labels (str)
obs[cluster_key + "_label"]     human-readable labels (if annotation_map provided)
uns["moranI"]                   DataFrame: genes × (I, pval_norm, pval_norm_fdr_bh)
uns["omicsage_spatial_cluster"] provenance dict
```

---

### 2.5 `spatial_deconvolve.py`

**Location:** `pipeline/modules/spatial/spatial_deconvolve.py`
**Phase 7 Session:** 4

**Purpose:** Deconvolve Visium spots into cell type abundances. Two methods
available: NNLS (default, fast, ~200 MB RAM) and cell2location (opt-in, deep
generative model, GPU recommended).

**Public API:**

```python
spatial_deconvolve(
    adata: AnnData,
    ref_adata: Optional[AnnData] = None,
    method: str = "nnls",                # "nnls" | "cell2location" | "none"
    cell_type_key: str = "cell_type_original",
    layer_ref: str = "counts",
    # NNLS-specific
    n_jobs: int = 4,
    target_sum: float = 10000,
    # cell2location-specific (ignored when method=nnls)
    batch_key_ref: Optional[str] = "donor_id",
    batch_key_st: Optional[str] = "patient",
    covariate_keys: Optional[list] = None,
    N_cells_per_location: int = 8,
    detection_alpha: int = 20,
    max_epochs_ref: int = 250,
    max_epochs_st: int = 30000,
    ...
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**Method selection:**
- `"nnls"` — default; non-negative least squares; fast, CPU-only, no GPU needed
- `"cell2location"` — opt-in; higher accuracy; requires GPU for large datasets
- `"none"` — explicit skip; appropriate for Xenium (cell-level data)
- `ref_adata=None` — always skips regardless of method

**Output keys added:**
```
obsm["q05_cell_abundance_w_sf"]      cell type abundances (n_spots × n_cell_types)
obsm["cell_type_abundances"]         canonical alias (same array)
obs[cell_type]                       per-cell-type abundance (one column per type)
uns["omicsage_spatial_deconvolve"]   provenance dict
```

---

### 2.6 `spatial_downstream.py`

**Location:** `pipeline/modules/spatial/spatial_downstream.py`
**Phase 7 Session:** 5

**Purpose:** All spatial downstream analyses in a single module — eight analyses
covering tissue niche identification, cell-type-resolved gene expression,
spatially variable gene characterisation, spatial interaction statistics,
ligand-receptor communication, and pathway enrichment.

*(Full API and analysis table unchanged — see previous version. Key parameters
and graceful skip conditions are the same.)*

**Provenance keys:** all `run_*` parameters, `timestamp`, `analyses` dict.

---

### 2.7 `spatial_impute.py`

**Location:** `pipeline/modules/spatial/spatial_impute.py`
**Phase 7 Session B** — *new*

**Purpose:** Impute full-transcriptome expression onto Visium/Visium HD spots
using a paired scRNA-seq reference. Fills the "missing genes" chapter: genes
not in the spatial panel are predicted from co-expression patterns in the
reference. Skips gracefully when no reference is configured.

**Public API:**

```python
spatial_impute(
    adata_spatial: AnnData,      # post-cluster checkpoint
    adata_sc: Optional[AnnData] = None,  # paired scRNA-seq reference
    method: str = "tangram",     # "tangram" | "gimvi"
    cell_type_key: str = "cell_type",
    n_top_genes: int = 2000,
    device: str = "cpu",
    tangram_mode: str = "clusters",   # "clusters" (memory-safe) | "cells" (GPU)
    max_cells_per_type: int = 500,    # subsample guard for cells mode only
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**Methods:**

| Method | Mode | RAM | Per-spot score | Notes |
|--------|------|-----|----------------|-------|
| `tangram` | `clusters` | ~100 MB | ✗ (not available) | Default; maps cell-type signatures |
| `tangram` | `cells` | ~2 GB+ | ✅ | Requires subsampling for large refs |
| `gimvi` | — | High, GPU rec | ✗ | Deep generative model via scvi-tools |

**Memory note — clusters mode:** In clusters mode, Tangram maps per-cell-type
mean signatures onto spots rather than individual cells. The mapping tensor is
`(n_cell_types × n_spots)` instead of `(n_cells × n_spots)` — memory is
negligible regardless of reference size. This is the recommended default for
production runs.

**Memory note — cells mode:** With 41k cells and 3k spots, the mapping tensor
is ~500 M float32 = ~2 GB before gradients. OOM-prone on WSL2. Use
`max_cells_per_type` to subsample (default 500 per type) if cells mode is needed.

**Tangram API notes** (verified from source):
- `tg.pp_adatas()` **mutates** both input objects — the module always works on copies
- `tg.map_cells_to_space(..., mode="clusters", cluster_label=cell_type_key)`
- `tg.project_genes(adata_map, adata_sc, cluster_label=cell_type_key)` — must pass
  `cluster_label` in clusters mode; the function collapses `adata_sc` to cluster-level
  internally before the obs-index check
- `project_genes` **returns** a new `AnnData` (spots × genes) — it does not
  mutate `ad_map`; the return value must be captured

**Storage contract:**
```
obsm["imputed_expression"]             float32 numpy array (n_spots × n_genes)
                                       h5py-safe; no DataFrame in obsm
obs["tangram_mapping_score"]           per-spot mapping score (cells mode only)
uns["omicsage_spatial_impute"]         provenance dict including "genes_imputed" list
```

Gene names are stored in `uns["omicsage_spatial_impute"]["outputs"]["genes_imputed"]`
and used by the report to reconstruct the DataFrame for visualization.

**Skip behaviour:**
- `adata_sc=None` → `skipped=True`, passthrough checkpoint written, report
  renders a skip notice
- `enabled: false` in config → same passthrough behaviour
- `sc_reference_path` empty/null in config → same

**Provenance keys:** `module`, `timestamp`, `method`, `skipped`, `skip_reason`,
`outputs.n_genes_imputed`, `outputs.n_spots`, `outputs.mean_mapping_score`,
`outputs.n_poor_spots`, `outputs.genes_imputed`, `outputs.cell_type_key`,
`outputs.device`, `outputs.tangram_mode`

---

## 3. Report Templates

All six reports follow the same structure:
- `_render_page(title, sections, timestamp)` → full HTML shell
- Content composed of `<section>` blocks inside `<main>`
- `.stat-grid`, `.stat-card`, `.fig-grid`, `.fig-wrap` CSS classes
- All figures base64-encoded PNG strings — fully self-contained
- Reports embeddable in the combined report via `<main>` extraction

---

### 3.1–3.5 (unchanged)

`spatial_qc_report.py`, `spatial_reduce_report.py`, `spatial_cluster_report.py`,
`spatial_deconvolve_report.py`, `spatial_downstream_report.py` — all unchanged
from Phase 7 Sessions 1–5. See previous version for section-level detail.

---

### 3.6 `spatial_impute_report.py`

**Location:** `reports/templates/spatial/spatial_impute_report.py`
**Phase 7 Session B** — *new*

**Public API:**
```python
generate_spatial_impute_report(
    adata: AnnData,           # must have uns["omicsage_spatial_impute"]
    output_path: str,
    dataset_id: str = "spatial",
    sc_ref_label: str = "",   # filename shown in Run Summary
) -> str                      # absolute path to written HTML
```

**Sections:**
1. **Run Summary** — stat cards: n_genes_imputed, method, mean mapping score
   (N/A in clusters mode — this is correct), n_spots, SC reference filename;
   explanatory note on mapping score interpretation
2. **Mapping Score Distribution** — histogram of per-spot Tangram mapping scores
   with red threshold line at 0.1 (cells mode only). In clusters mode, renders
   an explanatory paragraph instead of a skip note, explaining why scores are
   not available and how to get them
3. **Top Imputed Genes on Tissue** — `sq.pl.spatial_scatter` for top 5 genes by
   variance in imputed expression; squidpy-optional (skips gracefully if absent)
4. **Imputation Validation** — scatter of mean measured (log-normalised) vs mean
   imputed expression across up to 50 shared genes; Spearman r with colour-coded
   quality note (green ≥ 0.7, amber ≥ 0.4, red < 0.4). **Measured values are
   log-normalised** (`log1p(counts/lib_size × 10000)`) before computing Spearman
   to match the normalised scale of Tangram's imputed output

**Report DataFrame reconstruction:** The report reads
`obsm["imputed_expression"]` (numpy array) and
`uns["omicsage_spatial_impute"]["outputs"]["genes_imputed"]` (gene name list)
and reconstructs a DataFrame internally. Handles both live-run DataFrames
(if obsm still holds a DataFrame in memory) and checkpoint-reloaded numpy arrays.

---

## 4. Combined Report

**Location:** `reports/spatial_combined_report.py`

**Tab registry** (updated Session B):

| Filename | Tab label | Icon |
|----------|-----------|------|
| `spatial_qc_report.html` | QC | 🔬 |
| `spatial_reduce_report.html` | Reduce | 🔭 |
| `spatial_cluster_report.html` | Cluster | 🫧 |
| `spatial_deconvolve_report.html` | Deconvolve | 🧬 |
| `spatial_downstream_report.html` | Downstream | 🔗 |
| `spatial_impute_report.html` | Impute | 🧩 |

Progress bar now shows n_done / **6** steps complete.

All other combined report behaviour unchanged (tab extraction, keyboard nav,
lightbox, partial tab support).

---

## 5. Pipeline Runner

**Location:** `run_spatial_pipeline.py`

**Step order:** `ingest → qc → reduce → cluster → deconvolve → downstream → impute`

**Checkpoints:**

| Step | Checkpoint file |
|------|----------------|
| ingest | `01_ingested.h5ad` |
| qc | `02_qc.h5ad` |
| reduce | `03_reduced.h5ad` |
| cluster | `04_clustered.h5ad` |
| deconvolve | `05_deconvolved.h5ad` |
| downstream | `06_downstream.h5ad` |
| impute | `07_imputed.h5ad` |

**Impute step predecessor:** `"cluster"` (does not require deconvolution).
The runner loads the cluster checkpoint, runs imputation, generates the report,
and writes the imputed checkpoint. Writes a passthrough checkpoint (cluster
checkpoint copied + skip provenance) when imputation is disabled or no reference
is configured.

**Impute skip conditions** (checked in runner, not module):
- `config.spatial.impute.enabled: false`
- `config.spatial.impute.sc_reference_path` is null or empty string

**Runner reads from config:**
```yaml
spatial:
  impute:
    enabled: true
    method: tangram
    sc_reference_path: <path>
    n_top_genes: 2000
    cell_type_key: cell_type
    device: cpu
    tangram_mode: clusters      # "clusters" | "cells"
    max_cells_per_type: 500     # cells mode subsampling guard
```

**Usage (unchanged pattern):**
```bash
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --step impute
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step impute --force
```

---

## 6. Config Files

### `kuppe_heart.yaml` (updated Session B)

**Location:** `config/runs/kuppe_heart.yaml`

Added `spatial.impute` block:

```yaml
spatial:
  impute:
    enabled: true
    method: tangram
    sc_reference_path: data/benchmark/kuppe_snRNA_human_heart_2022_control.h5ad
    n_top_genes: 2000
    cell_type_key: "cell_type_original"   # matches deconvolve.cell_type_key
    device: cpu
    tangram_mode: clusters
    max_cells_per_type: 500
```

### `visium_hd_mouse_brain.yaml` (new — Session B)

**Location:** `config/runs/visium_hd_mouse_brain.yaml`

Configuration for 10x Visium HD Mouse Brain public dataset (used in the
Seurat Visium HD vignette). Key differences from Kuppe:

```yaml
spatial:
  source: data/benchmark/visium_hd_mouse_brain   # Space Ranger output dir
  spatial_type: visium_hd
  ingest:
    bin_size: 8             # 8µm bins (10x recommendation)
  qc:
    min_counts: 50          # lower threshold for smaller HD bins
    min_genes: 50
    mt_prefix: "mt-"        # mouse uses lowercase
```

**Download:** `https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain`
(free, no account required — direct download via `cf.10xgenomics.com`).

### `xenium_breast.yaml` (new — Session B)

**Location:** `config/runs/xenium_breast.yaml`

Configuration for 10x Xenium Human Breast Cancer public dataset. Key
differences from Visium:

```yaml
spatial:
  source: data/benchmark/xenium_breast
  spatial_type: xenium
  load_images: false          # morphology images not loaded by default
  qc:
    min_counts: 5             # panel-limited UMI
    min_genes: 5
    max_genes: 500            # targeted panel ceiling
  reduce:
    n_top_genes: 400          # <= panel size
    n_comps: 30               # fewer effective dimensions
  deconvolve:
    method: none              # cell-level data; deconvolution not applicable
  impute:
    enabled: false            # targeted panel; imputation less meaningful
```

**Download:** `https://www.10xgenomics.com/datasets/xenium-prime-5k-human-breast-cancer-ffpe`
(free, account registration required).

---

## 7. Tests

### `test_spatial_ingest_qc.py` (updated Session B)
**Location:** `tests/test_spatial_ingest_qc.py`

New test groups added in Session B:
- Auto-detection: `visium_hd` detected from `binned_outputs/` directory
- Auto-detection: `xenium` detected from `transcripts.parquet`
- Auto-detection priority: `visium_hd` wins over `visium` when both `spatial/`
  and `binned_outputs/` present
- Auto-detection priority: `xenium` wins over `visium` when both `spatial/`
  and `transcripts.parquet` present

For full Session B ingest format tests see `test_spatial_ingest_formats.py`.

---

### `test_spatial_ingest_formats.py` (new — Session B)

**Location:** `tests/test_spatial_ingest_formats.py`
**Covers:** Visium HD and Xenium loaders in `spatial_ingest.py`
**Strategy:** All tests mock `spatialdata-io` via `patch.dict(sys.modules, ...)` —
no real spatialdata-io install required in CI.

Key test groups:
- **Auto-detection:** all four new fingerprint tests (visium_hd, xenium, priority)
- **Visium HD loader:** output contract (obsm/uns/layers), bin_size 8 and 16 routing,
  missing table key → `KeyError` with message, `ImportError` without spatialdata-io
- **Xenium loader:** output contract, `metadata["platform"] == "xenium"`, missing
  table key → `KeyError`, `ImportError` without spatialdata-io
- **Registry:** `visium_hd` and `xenium` now show `"implemented"`;
  `merfish` and `codex` still `"planned"`
- **Still-planned formats:** `NotImplementedError` for merfish and codex
- **Public API routing:** `spatial_ingest(..., spatial_type="visium_hd")` and
  `spatial_ingest(..., spatial_type="xenium")` route correctly end-to-end

Total: 27 tests.

---

### `test_spatial_reduce.py`, `test_spatial_cluster.py`, `test_spatial_deconvolve.py`, `test_spatial_downstream.py`

Unchanged from Phase 7 Sessions 2–5. See previous version for detail.

---

### `test_spatial_impute.py` (new — Session B)

**Location:** `tests/test_spatial_impute.py`
**Covers:** `spatial_impute.py` + `spatial_impute_report.py`
**Strategy:** All Tangram calls mocked via `patch.object(m, "_TANGRAM_AVAILABLE", True)`
+ `patch.dict(sys.modules, {"tangram": mock_tg})`. No real Tangram install
required in CI. gimVI tested only for its `ImportError` guard.

Key test groups:
- **Skip path:** `adata_sc=None` skips cleanly; `skipped=True` in provenance;
  `inplace=False/True` semantics correct
- **Tangram path:** output contract (`obsm["imputed_expression"]` is `np.ndarray`,
  dtype `float32`); shape (spots × n_genes); provenance keys; `inplace=False`
  does not mutate input; invalid method → `ValueError`
- **gimVI ImportError guard:** raises `ImportError` when `_SCVI_AVAILABLE=False`
- **Tangram ImportError guard:** raises `ImportError` when `_TANGRAM_AVAILABLE=False`
- **Helpers:** `_ensure_counts_in_X` prefers `layers["counts"]`; `_select_overlap_hvgs`
  returns list within spatial genes; no-overlap → `ValueError`
- **Report:** renders without error (squidpy mocked out); skipped state renders
  skip notice; stat cards present; validation section present; numpy array obsm
  reconstructed from `uns` gene names correctly

Total: 20 tests.

---

## 8. Benchmark Data Download

**Location:** `scripts/download_spatial_benchmark.py`
**Phase 7 Session B** — *new*

**Purpose:** Downloads all spatial benchmark datasets needed for OmicSage testing
and validation. Kuppe datasets download automatically; 10x datasets print manual
instructions and validate the directory structure once placed.

**Usage:**
```bash
conda activate omicsage
python scripts/download_spatial_benchmark.py                        # download all
python scripts/download_spatial_benchmark.py --skip kuppe_visium    # skip one
python scripts/download_spatial_benchmark.py --validate-only        # check what's present
```

**Dataset registry:**

| Key | Dataset | Size | Method |
|-----|---------|------|--------|
| `kuppe_visium` | Kuppe et al. 2022 Visium (human heart) | ~1 GB | Auto (Figshare) |
| `kuppe_snrna` | Kuppe et al. 2022 snRNA-seq reference | ~2 GB | Auto (Figshare) |
| `visium_hd_mouse_brain` | 10x Visium HD Mouse Brain | ~3 GB | Manual (10x website) |
| `xenium_breast` | 10x Xenium Human Breast Cancer | ~8 GB | Manual (10x website, account required) |

**Validation checks:**
- Kuppe: file exists and `> 1 MB`
- Visium HD: directory contains `binned_outputs/` and `spatial/`
- Xenium: directory contains `experiment.xenium`, `cell_feature_matrix.h5`,
  `transcripts.parquet`

**Exit code:** `0` if all datasets ready or skipped; `1` if any manual download
is required or a download failed (CI-compatible).

---

## 9. AnnData State Reference

The table below shows the complete state of `adata` after each pipeline step.

| Key | ingest | QC | reduce | cluster | deconvolve | downstream | impute |
|-----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `X` | raw | raw | norm+log1p | norm+log1p | norm+log1p | norm+log1p | norm+log1p |
| `layers["counts"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obsm["spatial"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["spatial"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["total_counts"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["qc_pass"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `var["highly_variable"]` | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obsm["X_pca"]` | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obsp["spatial_connectivities"]` | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["spatial_cluster"]` | — | — | — | ✅ | ✅ | ✅ | ✅ |
| `uns["moranI"]` | — | — | — | ✅¹ | ✅ | ✅ | ✅ |
| `obsm["q05_cell_abundance_w_sf"]` | — | — | — | — | ✅² | ✅² | ✅² |
| `obs[cell_type columns]` | — | — | — | — | ✅² | ✅² | ✅² |
| `obs["region_cluster"]` | — | — | — | — | — | ✅² | ✅² |
| `obsm["X_umap_celltype"]` | — | — | — | — | — | ✅² | ✅² |
| `uns["celltype_marker_genes"]` | — | — | — | — | — | ✅² | ✅² |
| `uns["*_nhood_enrichment"]` | — | — | — | — | — | ✅² | ✅² |
| `uns["*_co_occurrence"]` | — | — | — | — | — | ✅² | ✅² |
| `uns["*_ligrec"]` | — | — | — | — | — | ✅² | ✅² |
| `uns["svg_gsea"]` | — | — | — | — | — | ✅ | ✅ |
| `obsm["imputed_expression"]` | — | — | — | — | — | — | ✅³ |
| `obs["tangram_mapping_score"]` | — | — | — | — | — | — | ✅⁴ |
| `uns["omicsage_spatial_ingest"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_qc"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_reduce"]` | — | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_cluster"]` | — | — | — | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_deconvolve"]` | — | — | — | — | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_downstream"]` | — | — | — | — | — | ✅ | ✅ |
| `uns["omicsage_spatial_impute"]` | — | — | — | — | — | — | ✅ |

¹ Only if `run_svg=True`
² Only if `ref_adata` was provided to `spatial_deconvolve`
³ float32 numpy array; `uns["omicsage_spatial_impute"]["outputs"]["genes_imputed"]`
  holds column names
⁴ Only in `tangram_mode="cells"`; not available in clusters mode

`*` = `dominant_celltype_key` prefix (default: `dominant_cell_type`)

---

## 10. Key Decisions

| Decision | Choice | Reason |
|----------|--------|--------|
| Normalization | `normalize_total` + `log1p` | SCTransform is R-only; Python/scanpy standard |
| HVG selection | `flavor="seurat"`, `n_top_genes=3000` | Safe for RNA; 3000 standard for Visium |
| Spatial graph | `coord_type=None` (auto-detect) | squidpy auto-selects grid; safer than hardcoding |
| Clustering graph | KNN from PCA, NOT spatial graph | Cluster by transcriptomic similarity; spatial graph reserved for SVGs |
| SVG p-values | `n_perms=None` (analytical) | Fast; permutation testing available but slow |
| Deconvolution default | NNLS | Memory-safe; cell2location opt-in |
| Xenium deconvolution | `method: none` | Cell-level data; deconvolution not applicable |
| Visium HD bin default | 8µm | 10x recommendation; 2µm is too sparse |
| Visium HD / Xenium loader | `spatialdata-io` | Both `sq.read.visium_hd` and `sq.read.xenium` do not exist in squidpy; spatialdata-io is the correct package |
| Visium HD / Xenium SpatialData→AnnData | Extract `sdata.tables[key]` | Tables already have `obsm["spatial"]` set — no extra extraction needed |
| uns["spatial"] stub for new formats | Minimal stub built in loader | Downstream squidpy tools validate for `uns["spatial"]`; imaging data not loaded by default |
| Imputed expression storage | `float32` numpy array in `obsm` | AnnData rejects plain Python strings in obsm; DataFrames fail h5py due to mixed-type string columns |
| Gene names for imputed expression | Stored in `uns` provenance | Decouples gene identity from the array; compatible with h5ad checkpoint round-trip |
| Tangram mode default | `"clusters"` | `"cells"` mode OOM on 41k cell reference (41k × 3k spots tensor); clusters mode uses n_cell_types × n_spots — negligible |
| Tangram `project_genes` | Must pass `cluster_label` + capture return value | In clusters mode: `cluster_label` required for obs-index alignment; `project_genes` returns new AnnData, does not mutate `ad_map` |
| Validation scatter normalization | log-normalise measured before Spearman | Imputed values are on normalised scale; raw counts vs normalised gives artificially low r |
| Mapping scores in clusters mode | Not computed | `ad_map.obs` rows are cell types, not spots; `tg_score` cannot map back to spots |
| Imputation for Xenium | `enabled: false` default | Xenium is already a targeted panel; imputing genes not in the panel has less value than for Visium |
| Gene ID consistency | `_select_overlap_hvgs` uses `set` intersection of `var_names` | Ensures no mismatched symbols vs ENSEMBL IDs; raises `ValueError` with clear message if overlap is empty |
| Annotation map | Optional, no hardcoded maps | Cluster numbering is non-deterministic across runs |
| `inplace` default | `False` | Safe default |
| Report structure | `<header>/<main>/<footer>` with `<section>` blocks | Required for combined report tab extraction |
| GSEA dtype sanitisation | `_sanitize_gsea_df()` before `uns["svg_gsea"]` | gseapy returns numeric columns as `object` dtype; `write_h5ad` crashes |
| Downstream predecessor | `"cluster"` (minimum); upgrades to deconvolve if available | Allows running on any dataset without scRNA-seq reference |
| Impute predecessor | `"cluster"` (does not require deconvolution) | Imputation is independent of cell type deconvolution |
