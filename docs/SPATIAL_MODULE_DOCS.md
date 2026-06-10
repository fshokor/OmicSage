# OmicSage — Spatial Transcriptomics Module Documentation
> Phase 7 — Spatial Transcriptomics (Sessions 1–5)
> Last updated: June 2026

This document covers every script produced in Phase 7: pipeline modules,
report generators, the combined report, the runner, the config, and the tests.

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
3. [Report Templates](#3-report-templates)
   - [spatial_qc_report.py](#31-spatial_qc_reportpy)
   - [spatial_reduce_report.py](#32-spatial_reduce_reportpy)
   - [spatial_cluster_report.py](#33-spatial_cluster_reportpy)
   - [spatial_deconvolve_report.py](#34-spatial_deconvolve_reportpy)
   - [spatial_downstream_report.py](#35-spatial_downstream_reportpy)
4. [Combined Report](#4-combined-report)
5. [Pipeline Runner](#5-pipeline-runner)
6. [Config File](#6-config-file)
7. [Tests](#7-tests)
8. [AnnData State Reference](#8-anndata-state-reference)
9. [Key Decisions](#9-key-decisions)

---

## 1. Pipeline Architecture

```
spatial_ingest.py
      │  standard AnnData: obsm["spatial"], uns["spatial"], layers["counts"]
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
**Phase 7 Session:** 1 (updated Session 4)

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
    load_images: bool = True,
    inplace: bool = False,
) -> tuple[AnnData, dict]

list_supported_types() -> dict[str, str]  # maps type → "implemented"|"planned"
```

**Supported technologies:**

| Type | Status | Source format |
|------|--------|---------------|
| `visium` | ✅ Implemented | Space Ranger output directory |
| `h5ad` | ✅ Implemented | Pre-built `.h5ad` file |
| `benchmark` | ✅ Implemented | `squidpy.datasets.visium_hne_adata()` |
| `visium_hd` | 🔲 Planned | `binned_outputs/` directory |
| `xenium` | 🔲 Planned | `transcripts.parquet` |
| `merfish` | 🔲 Planned | `cell_by_gene.csv` |
| `codex` | 🔲 Planned | `.csv` protein markers |

**Auto-detection fingerprints** (checked in order):
1. `source == "benchmark"` → benchmark
2. `source.endswith(".h5ad")` → h5ad
3. Directory contains `transcripts.parquet` → xenium
4. Directory contains `cell_by_gene.csv` → merfish
5. Directory contains `binned_outputs/` → visium_hd
6. Directory contains `spatial/` subfolder → visium
7. File ends with `.csv` → codex

**h5ad loading contract** (updated Session 4 for Kuppe data):
- If `layers["counts"]` absent → saves `X.copy()` to `layers["counts"]`
- If `var["gene_ids"]` present → swaps `var_names` to ENSEMBL IDs, preserves symbols in `var["feature_name"]`
- If any `MT-` genes present → strips them into `obsm["MT"]` and removes from `X`

**Output AnnData contract** (all implemented types):
```
obsm["spatial"]                             coordinates (n_obs, 2)
uns["spatial"][library_id]["images"]        tissue images
uns["spatial"][library_id]["scalefactors"]  scale factors
layers["counts"]                            raw integer counts
uns["omicsage_spatial_ingest"]              provenance dict
```

**Provenance keys:** `source`, `spatial_type`, `technology_notes`, `counts_file`, `library_id`, `load_images`, `n_obs`, `n_vars`, `timestamp`

---

### 2.2 `spatial_qc.py`

**Location:** `pipeline/modules/spatial/spatial_qc.py`
**Phase 7 Session:** 1

**Purpose:** Compute QC metrics for Visium spots and optionally filter
low-quality spots based on count/gene/MT% thresholds.

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
total_counts          total UMI counts per spot
n_genes_by_counts     number of detected genes per spot
pct_counts_mt         mitochondrial gene percentage
qc_pass               bool — True if spot passes all thresholds
```

**Provenance keys:** all threshold parameters, `n_spots_before`, `n_spots_after`, `n_spots_removed`, per-filter removal counts, summary stats per metric.

**Note on Kuppe Visium data:** After `spatial_ingest` strips MT genes into `obsm["MT"]`, `pct_counts_mt` will be 0 for all spots. The `max_mt_pct` threshold is vacuously satisfied — this is correct behaviour.

---

### 2.3 `spatial_reduce.py`

**Location:** `pipeline/modules/spatial/spatial_reduce.py`
**Phase 7 Session:** 2

**Purpose:** Normalize counts, select highly variable genes, compute PCA,
and build the spatial neighbours graph. This is the bridge between raw
QC-filtered counts and Leiden clustering.

**Public API:**

```python
spatial_reduce(
    adata: AnnData,
    n_top_genes: int = 3000,
    n_comps: int = 50,
    n_neighbors: int = 6,                # 6 = Visium hex grid
    coord_type: Optional[str] = None,    # None = auto-detect grid
    normalize_total: bool = True,
    target_sum: float = 1e4,
    log1p: bool = True,
    flavor: str = "seurat",
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**What it does:**
1. **Raw count restoration** — if `layers["counts"]` present, resets `X` to raw before normalizing
2. **Benchmark detection** — if `spatial_type == "benchmark"` and no `layers["counts"]`, skips normalization (data already processed)
3. **Normalization** — `sc.pp.normalize_total` → `sc.pp.log1p` (skipped if `normalize_total=False`)
4. **HVG selection** — `sc.pp.highly_variable_genes(flavor="seurat", n_top_genes=3000)`
5. **PCA** — `sc.pp.pca(use_highly_variable=True, svd_solver="arpack")`; `n_comps` capped at `min(n_comps, n_hvgs-1, n_obs-1)`
6. **Spatial neighbours** — `sq.gr.spatial_neighbors(coord_type=None)` — auto-detects grid when `uns["spatial"]` present

**Key decision — `coord_type=None`:** squidpy auto-selects `"grid"` when `uns["spatial"]` is present (standard Visium). Safer than hardcoding `"grid"` as it will auto-adapt for future technology types.

**Output keys added:**
```
X                              normalized + log1p counts
layers["counts"]               raw counts (preserved)
var["highly_variable"]         HVG boolean mask
var["means"], var["dispersions_norm"]
obsm["X_pca"]                  PCA embedding (n_obs × n_comps)
uns["pca"]["variance_ratio"]   explained variance per PC
obsp["spatial_connectivities"] spatial adjacency matrix
obsp["spatial_distances"]      spatial distances
uns["spatial_neighbors"]       spatial graph metadata
uns["omicsage_spatial_reduce"] provenance dict
```

**Provenance keys:** all parameters, `skipped_normalization`, `n_hvgs`, `n_comps_computed`, `pca_variance_ratio_top10`, `pca_cumulative_variance_top10`, `spatial_graph_n_edges`, `spatial_graph_mean_neighbors`.

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
    n_neighbors: int = 15,               # KNN graph for Leiden (NOT spatial graph)
    n_pcs: int = 30,
    random_state: int = 0,
    cluster_key: str = "spatial_cluster",
    annotation_map: Optional[dict] = None,  # cluster_id → label
    run_svg: bool = True,
    svg_n_genes: Optional[int] = None,   # None = all HVGs
    svg_n_jobs: int = 1,
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**What it does:**
1. **KNN graph** — `sc.pp.neighbors(use_rep="X_pca")` — transcriptomic similarity graph (NOT the spatial graph)
2. **Leiden clustering** — `sc.tl.leiden(resolution=...)` → `obs[cluster_key]`
3. **Optional annotation** — maps cluster IDs to human-readable labels via `annotation_map`; writes `obs[cluster_key + "_label"]`
4. **Moran's I SVGs** — `sq.gr.spatial_autocorr(mode="moran", n_perms=None)` on HVGs only; uses analytical p-values for speed; results sorted descending by I score

**Key design decision:** Leiden uses the transcriptomic KNN graph, NOT the spatial graph. The spatial graph (`obsp["spatial_connectivities"]`) is reserved for Moran's I. This mirrors sc-best-practices: cluster by gene expression similarity, detect spatial patterns separately.

**Output keys added:**
```
obs[cluster_key]               Leiden cluster labels (str)
obs[cluster_key + "_label"]    human-readable labels (if annotation_map provided)
uns["moranI"]                  DataFrame: genes × (I, pval_norm, pval_norm_fdr_bh)
                                sorted descending by I score
uns["omicsage_spatial_cluster"] provenance dict
```

**Provenance keys:** all parameters, `n_clusters`, `cluster_sizes` (dict), `n_annotated_spots`, `n_genes_tested`, `n_significant_fdr05`, `top5_svg`.

---

### 2.5 `spatial_deconvolve.py`

**Location:** `pipeline/modules/spatial/spatial_deconvolve.py`
**Phase 7 Session:** 4

**Purpose:** Deconvolve Visium spots into cell type abundances using
cell2location. Each ~55µm Visium spot contains multiple cells — deconvolution
estimates what proportion of each cell type is present in each spot.

**Requires:** A paired scRNA-seq reference AnnData from the same tissue type.
When `ref_adata=None`, the module skips gracefully and records `skipped=True`.

**Public API:**

```python
spatial_deconvolve(
    adata: AnnData,
    ref_adata: Optional[AnnData] = None,
    cell_type_key: str = "cell_type_original",  # obs col in ref
    batch_key_ref: Optional[str] = "donor_id",  # batch key in ref
    batch_key_st: Optional[str] = "patient",    # batch key in spatial
    covariate_keys: Optional[list] = None,      # e.g. ["assay"]
    layer_ref: str = "counts",                  # raw counts layer in ref
    N_cells_per_location: int = 8,              # expected cells per spot
    detection_alpha: int = 20,
    max_epochs_ref: int = 250,
    max_epochs_st: int = 30000,
    batch_size_ref: int = 2500,
    cell_count_cutoff: int = 5,
    cell_percentage_cutoff2: float = 0.03,
    nonz_mean_cutoff: float = 1.12,
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

All cell2location parameters are arguments — different tissues require
different values. Defaults are tuned for the Kuppe human heart dataset.

**What it does (when `ref_adata` provided):**
1. **Prepare reference** — resets `ref.X` to raw counts from `layers[layer_ref]`
2. **Gene selection** — `c2l.utils.filtering.filter_genes` with permissive thresholds
3. **Subset to shared genes** — intersects spatial and reference gene sets; raises if empty
4. **Reference model** — fits `RegressionModel` → learns cell type expression signatures (`inf_aver`)
5. **Spatial model** — fits `Cell2location` with `N_cells_per_location` and `detection_alpha`
6. **Export posterior** — 1000 posterior samples → `obsm["q05_cell_abundance_w_sf"]`
7. **Write per-cell-type columns** — `obs[cell_type]` for each cell type (5% quantile abundance)

**Device:** GPU if available, CPU fallback — handled automatically by PyTorch/scvi-tools. No explicit `use_gpu` parameter needed.

**Benchmark note:** The squidpy benchmark dataset (`visium_hne_adata`) has no paired scRNA-seq reference. Pass `ref_adata=None` for benchmark runs — deconvolution is skipped, all other pipeline steps complete normally.

**Output keys added:**
```
obsm["q05_cell_abundance_w_sf"]      cell type abundances (n_spots × n_cell_types)
obs[cell_type]                       per-cell-type 5% quantile abundance (one column per type)
uns["cell2location_mod"]             cell2location model metadata
uns["omicsage_spatial_deconvolve"]   provenance dict
```

**Provenance keys:** `skipped`, `skip_reason` (when skipped), all parameters, `n_cell_types`, `cell_type_names`, `n_shared_genes`, `n_genes_after_selection`, `n_spots`.

**Kuppe dataset specifics:**
- `cell_type_key = "cell_type_original"` (11 cell types)
- `batch_key_ref = "donor_id"`
- `batch_key_st = "patient"` (4 patients: P1, P7, P8, P17)
- `covariate_keys = ["assay"]`
- `N_cells_per_location = 8`

---

### 2.6 `spatial_downstream.py`

**Location:** `pipeline/modules/spatial/spatial_downstream.py`
**Phase 7 Session:** 5

**Purpose:** All spatial downstream analyses in a single module — seven analyses
covering tissue niche identification, cell-type-resolved gene expression,
spatially variable gene characterisation, spatial interaction statistics,
ligand-receptor communication, and pathway enrichment. Every analysis is
independently gated: if its required inputs are absent, it records
`skipped=True` in provenance and the pipeline continues without error.

**Public API:**

```python
spatial_downstream(
    adata: AnnData,
    # Region clustering (sc-best-practices §32.3.4.2)
    run_region_clustering: bool = True,
    region_resolution: float = 0.5,
    region_n_neighbors: int = 15,
    # Cell-type specific gene expression
    run_celltype_expression: bool = True,
    n_marker_genes: int = 20,
    # Cell-type specific SVGs
    run_celltype_svg: bool = True,
    svg_n_genes: Optional[int] = None,    # None = all HVGs
    # Co-occurrence
    run_co_occurrence: bool = True,
    co_occurrence_interval: Optional[list] = None,
    # Neighbourhood enrichment
    run_nhood_enrichment: bool = True,
    n_perms_nhood: int = 1000,
    # Ligand-receptor
    run_ligrec: bool = True,
    ligrec_n_perms: int = 1000,
    ligrec_organism: str = "human",       # "human" or "mouse"
    # SVG pathway enrichment
    run_svg_gsea: bool = True,
    svg_gsea_gene_sets: str = "GO_Biological_Process_2023",
    svg_gsea_organism: str = "Human",     # "Human" or "Mouse"
    # Shared
    dominant_celltype_key: str = "dominant_cell_type",
    n_jobs: int = 1,
    inplace: bool = False,
) -> tuple[AnnData, dict]
```

**Analyses — Tier 1** (require only `uns["moranI"]`, available after `spatial_cluster`):

| Analysis | Method | Output key |
|----------|--------|-----------|
| SVG pathway enrichment | `gseapy.prerank` on Moran's I ranking | `uns["svg_gsea"]` |

**Analyses — Tier 2** (require deconvolution outputs):

| Analysis | Method | Output key |
|----------|--------|-----------|
| Region clustering | Leiden on `obsm["q05_cell_abundance_w_sf"]` KNN graph | `obs["region_cluster"]`, `obsm["X_umap_celltype"]` |
| Cell-type expression | Vectorised rank-based Spearman, all spots × all genes | `uns["celltype_marker_genes"]` |
| Cell-type specific SVGs | Moran's I on above-median abundance spot subsets | `uns["celltype_svg"]` |
| Co-occurrence | `sq.gr.co_occurrence` across distance intervals | `uns["{key}_co_occurrence"]` |
| Neighbourhood enrichment | `sq.gr.nhood_enrichment`, 1000 permutations | `uns["{key}_nhood_enrichment"]` |
| Ligand-receptor | `sq.gr.ligrec` via OmniPath database | `uns["{key}_ligrec"]` |

`{key}` = `dominant_celltype_key` (default: `"dominant_cell_type"`), following squidpy's naming convention.

**Graceful skip conditions:**

| Analysis | Skips when |
|----------|-----------|
| Region clustering | `obsm["q05_cell_abundance_w_sf"]` absent |
| Cell-type expression | `obsm["q05_cell_abundance_w_sf"]` absent, or no cell types resolved |
| Cell-type SVGs | `uns["moranI"]` or `obsm["q05_cell_abundance_w_sf"]` absent |
| Co-occurrence | `obs["dominant_cell_type"]` or `obsm["spatial"]` absent |
| Neighbourhood enrichment | `obs["dominant_cell_type"]` absent |
| Ligand-receptor | `obs["dominant_cell_type"]` absent, squidpy not installed, or `sq.gr.ligrec` unavailable |
| SVG GSEA | `uns["moranI"]` absent, or gseapy not installed |

**What it does — selected implementation notes:**

1. **Region clustering** — builds a KNN graph in cell-type abundance space (`sc.pp.neighbors` with `key_added="neighbors_celltype"`) so it does not overwrite the gene-expression KNN graph from earlier steps. Leiden and UMAP run on that graph. UMAP result stored at `obsm["X_umap_celltype"]`; any pre-existing `obsm["X_umap"]` is saved and restored.

2. **Cell-type expression** — vectorised rank-based Spearman computed as a matrix operation rather than per-gene scipy calls. Gene names are mapped from ENSEMBL IDs to symbols via `var["feature_name"]` before storage.

3. **Cell-type SVGs** — each cell-type subset recomputes `sq.gr.spatial_neighbors(n_neighs=6)` because the original `obsp["spatial_connectivities"]` is invalid after subsetting. Uses `n_perms=None` (analytical p-values) for speed.

4. **Ligand-receptor** — makes a temporary `adata.copy()` with `var_names` remapped from ENSEMBL IDs to `var["feature_name"]` symbols before calling `sq.gr.ligrec`. The result is copied back; the original `var_names` are never modified. Requires internet access to query OmniPath at runtime.

5. **SVG GSEA** — `_sanitize_gsea_df()` casts the five known numeric columns (`ES`, `NES`, `NOM p-val`, `FDR q-val`, `FWER p-val`) to `float64` and all others to `str` before storing in `uns`. Without this, `adata.write_h5ad()` crashes because h5py cannot serialise Python floats as variable-length strings when columns have `object` dtype.

**Provenance keys:** all `run_*` parameters, `timestamp`, `analyses` dict — one entry per analysis with `skipped`, `reason` (when skipped), and analysis-specific stats (e.g. `n_regions`, `n_cell_types`, `n_pathways`, `n_significant`).

---

## 3. Report Templates

All five reports follow the same structure as the RNA pipeline reports:
- `_render_page(title, sections, timestamp)` → full HTML shell with `<header>`, `<main>`, `<footer>`
- Content inside `<main>` is composed of `<section>` blocks
- Sections use `.stat-grid`, `.stat-card`, `.fig-grid`, `.fig-wrap` CSS classes
- All figures are base64-encoded PNG strings embedded directly in the HTML
- Reports are self-contained single files — no external dependencies

**Why this pattern?** The `generate_spatial_combined_report()` function
extracts content from `<main>` tags to assemble the tabbed combined report.
All spatial reports must use `<main>` for this to work correctly.

---

### 3.1 `spatial_qc_report.py`

**Location:** `reports/templates/spatial/spatial_qc_report.py`

**Public API:**
```python
generate_spatial_qc_report(
    adata: AnnData,           # must have uns["omicsage_spatial_qc"]
    output_path: str,
    dataset_id: str = "spatial",
) -> str                      # absolute path to written HTML
```

**Sections:**
1. **Run Summary** — stat cards (spots in/kept/removed, pass rate, genes); QC metric summary table (mean/median/std/min/max); filter thresholds applied table
2. **QC Metric Distributions** — violin plots: total_counts, n_genes_by_counts, pct_counts_mt with threshold lines
3. **Spatial Distribution of QC Metrics** — `sq.pl.spatial_scatter` for total_counts and pct_counts_mt overlaid on tissue image
4. **Filter Breakdown** — bar chart of spots removed per filter criterion

---

### 3.2 `spatial_reduce_report.py`

**Location:** `reports/templates/spatial/spatial_reduce_report.py`

**Public API:**
```python
generate_spatial_reduce_report(
    adata: AnnData,           # must have uns["omicsage_spatial_reduce"]
    output_path: str,
    dataset_id: str = "spatial",
) -> str
```

**Sections:**
1. **Run Summary** — stat cards (spots, total genes, HVGs, PCA components, variance explained, spatial edges, mean neighbours); parameters table; normalization-skipped note if applicable
2. **Highly Variable Genes** — scatter: mean expression vs normalized dispersion, HVGs highlighted in orange
3. **PCA** — elbow plot (per-PC + cumulative variance); PCA scatter coloured by total_counts
4. **Spatial Neighbours Graph** — histogram of neighbours per spot (expected ~6 for Visium grid)

---

### 3.3 `spatial_cluster_report.py`

**Location:** `reports/templates/spatial/spatial_cluster_report.py`

**Public API:**
```python
generate_spatial_cluster_report(
    adata: AnnData,           # must have uns["omicsage_spatial_cluster"]
    output_path: str,
    dataset_id: str = "spatial",
) -> str
```

**Sections:**
1. **Run Summary** — stat cards (spots, clusters, resolution, SVGs if run); top 5 SVG gene badges; parameters table
2. **Clustering** — spatial scatter coloured by Leiden cluster; horizontal bar chart of spots per cluster
3. **Spatially Variable Genes (Moran's I)** — horizontal bar chart of top 20 SVG scores (red = FDR < 0.05); spatial scatter of top 3 SVG genes; table of top 20 SVGs with I score, p-value, FDR-adjusted p-value

The SVG section is omitted entirely when `run_svg=False`.

---

### 3.4 `spatial_deconvolve_report.py`

**Location:** `reports/templates/spatial/spatial_deconvolve_report.py`

**Public API:**
```python
generate_spatial_deconvolve_report(
    adata: AnnData,           # must have uns["omicsage_spatial_deconvolve"]
    output_path: str,
    dataset_id: str = "spatial",
) -> str
```

**Sections:**
1. **Run Summary** — stat cards (spots, cell types or "Skipped", shared genes); skipped notice with reason if applicable; cell type badges
2. **Cell Type Abundances** *(only when not skipped)* — mean abundance bar chart; dominant cell type per spot (spatial scatter); top 6 cell types spatial distribution (2×3 grid)
3. **Parameters Used** *(only when not skipped)* — full parameters table

When `skipped=True`, only the Run Summary section is shown with a yellow
warning note explaining that a paired scRNA-seq reference is required.

---

### 3.5 `spatial_downstream_report.py`

**Location:** `reports/templates/spatial/spatial_downstream_report.py`

**Public API:**
```python
generate_spatial_downstream_report(
    adata: AnnData,           # must have uns["omicsage_spatial_downstream"]
    output_path: str,
    dataset_id: str = "spatial",
    dominant_celltype_key: str = "dominant_cell_type",
) -> str
```

**Sections:**
1. **Run Summary** — stat cards (region clusters, cell types profiled, cell types SVG-tested, SVG pathways); analysis status table showing run/skipped for all 7 analyses
2. **Region Clustering** — spatial scatter coloured by `obs["region_cluster"]`; UMAP of cell type composition space (`obsm["X_umap_celltype"]`) if available
3. **Cell-type Marker Genes** — table: top 10 Spearman-correlated genes per cell type
4. **Cell-type Specific SVGs** — table: top 5 Moran's I SVGs per cell type subset
5. **Spatial Co-occurrence** — `sq.pl.co_occurrence` line plot for most abundant cell type
6. **Neighbourhood Enrichment** — `sq.pl.nhood_enrichment` z-score heatmap (`method="average"`)
7. **Ligand-Receptor Communication** — `sq.pl.ligrec` dotplot at `alpha=0.001`; significant interaction count note
8. **SVG Pathway Enrichment** — NES bar chart (top 10, red = positive / blue = negative); top 20 pathway table with FDR significance stars

Every section renders a "not run / data not available" note when the corresponding analysis was skipped, so the report is always complete regardless of which analyses ran.

---

## 4. Combined Report

**Location:** `reports/spatial_combined_report.py`

**Purpose:** Assembles all five step HTML reports into a single self-contained
tabbed HTML file. Identical tab UI to `reports/combined_report.py` (RNA pipeline).

**Public API:**
```python
generate_spatial_combined_report(
    reports_dir: Path,
    dataset_name: str = "OmicSage Spatial Analysis",
    output_path: Optional[Path] = None,   # default: reports_dir/00_spatial_combined_report.html
) -> str                                  # absolute path to written file
```

**Tab registry:**

| Filename | Tab label | Icon |
|----------|-----------|------|
| `spatial_qc_report.html` | QC | 🔬 |
| `spatial_reduce_report.html` | Reduce | 🔭 |
| `spatial_cluster_report.html` | Cluster | 🫧 |
| `spatial_deconvolve_report.html` | Deconvolve | 🧬 |
| `spatial_downstream_report.html` | Downstream | 🔗 |

**Behaviour:**
- Only tabs for reports that exist on disk are shown
- Progress bar shows n_done / 5 steps complete
- Keyboard navigation: left/right arrow keys switch tabs
- Supports both bare filenames and dataset-prefixed filenames (e.g. `kuppe_heart_spatial_qc_report.html`)

**Standalone CLI:**
```bash
python -m reports.spatial_combined_report \
  --reports-dir outputs/kuppe_heart \
  --dataset-name "Kuppe Heart 2022" \
  --output outputs/kuppe_heart/00_spatial_combined_report.html
```

---

## 5. Pipeline Runner

**Location:** `run_spatial_pipeline.py`

**Purpose:** Orchestrates all spatial pipeline steps with checkpointing,
step selection, and combined report generation. Matches the RNA
`run_pipeline.py` pattern exactly.

**Usage:**
```bash
# Run all steps
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml

# Run all steps (explicit)
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --step all

# Stop at a checkpoint (inclusive)
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --to-step cluster

# Resume from a checkpoint (inclusive)
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step reduce

# Run a specific range
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step qc --to-step cluster

# Run exactly one step
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --step cluster

# Force re-run even if checkpoint exists
python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step reduce --force
```

**Step order:** `ingest → qc → reduce → cluster → deconvolve → downstream`

**Checkpointing:** Every step writes its output to `output_dir/NN_<step>.h5ad`.
If the file already exists, the step is skipped (cached). Use `--force` to
override, or `--from-step` to force re-execution from a given step onward.

| Step | Checkpoint file |
|------|----------------|
| ingest | `01_ingested.h5ad` |
| qc | `02_qc.h5ad` |
| reduce | `03_reduced.h5ad` |
| cluster | `04_clustered.h5ad` |
| deconvolve | `05_deconvolved.h5ad` |
| downstream | `06_downstream.h5ad` |

**Downstream step predecessor logic:** `STEP_PREDECESSOR["downstream"] = "cluster"` (minimum required). The runner's `resolve_input` upgrades to `05_deconvolved.h5ad` automatically if it exists on disk, so running deconvolution first always gives the richer result without any extra CLI flags.

**Config lookup:** All step parameters are read from the YAML config under
`spatial.<step>`. CLI overrides are not supported (use `--force` + config edits).

**Combined report:** Generated automatically at the end of every run, covering
all step reports that exist in `output_dir`.

---

## 6. Config File

**Location:** `config/runs/kuppe_heart.yaml`
**Dataset:** Kuppe et al. 2022 human myocardial infarction (control samples)

**Top-level keys:**

| Key | Value | Description |
|-----|-------|-------------|
| `dataset_id` | `kuppe_heart` | Used in output paths and report titles |
| `dataset_name` | *(optional)* | Human-readable name for combined report header |
| `output_dir` | `outputs/kuppe_heart` | All checkpoints and reports written here |

**`spatial` section keys:**

| Key | Default | Description |
|-----|---------|-------------|
| `source` | *(required)* | Path to h5ad file or `"benchmark"` |
| `spatial_type` | `"auto"` | `"h5ad"`, `"visium"`, `"benchmark"`, etc. |
| `load_images` | `true` | Load tissue images from h5ad |

**`spatial.qc` keys:**

| Key | Default | Description |
|-----|---------|-------------|
| `min_counts` | `500` | Minimum total UMI per spot |
| `max_counts` | `100000` | Maximum total UMI per spot |
| `min_genes` | `200` | Minimum genes detected per spot |
| `max_genes` | `10000` | Maximum genes detected per spot |
| `max_mt_pct` | `20.0` | Maximum mitochondrial % |
| `mt_prefix` | `"MT-"` | `"MT-"` human, `"mt-"` mouse |
| `filter_spots` | `true` | Whether to remove failing spots |

**`spatial.reduce` keys:**

| Key | Default | Description |
|-----|---------|-------------|
| `n_top_genes` | `3000` | HVGs to select |
| `n_comps` | `50` | PCA components |
| `n_neighbors` | `6` | Spatial graph neighbours (6 = Visium hex grid) |
| `coord_type` | `null` | `null` = auto-detect; `"grid"` or `"generic"` to override |
| `normalize_total` | `true` | Apply normalize_total |
| `target_sum` | `10000` | Normalization target (CPM-like) |
| `log1p` | `true` | Apply log1p after normalization |
| `flavor` | `"seurat"` | HVG flavor: `"seurat"` or `"cell_ranger"` |

**`spatial.cluster` keys:**

| Key | Default | Description |
|-----|---------|-------------|
| `resolution` | `0.5` | Leiden resolution |
| `n_neighbors` | `15` | KNN graph neighbours (for Leiden, NOT spatial) |
| `n_pcs` | `30` | PCA components to use for KNN |
| `random_state` | `0` | Random seed |
| `run_svg` | `true` | Run Moran's I SVG detection |
| `svg_n_genes` | `null` | Genes to test (`null` = all HVGs) |
| `annotation_map` | `null` | `cluster_id: label` dict or `null` |

**`spatial.deconvolve` keys:**

| Key | Default | Description |
|-----|---------|-------------|
| `ref_path` | *(optional)* | Path to scRNA-seq reference h5ad |
| `cell_type_key` | `"cell_type_original"` | obs column in reference with cell types |
| `batch_key_ref` | `"donor_id"` | Batch key in reference |
| `batch_key_st` | `"patient"` | Batch key in spatial data |
| `covariate_keys` | `["assay"]` | Categorical covariates in reference |
| `layer_ref` | `"counts"` | Raw counts layer in reference |
| `N_cells_per_location` | `8` | Expected cells per Visium spot |
| `detection_alpha` | `20` | Technical variability parameter |
| `max_epochs_ref` | `250` | Reference model training epochs |
| `max_epochs_st` | `30000` | Spatial model training epochs |
| `batch_size_ref` | `2500` | Reference model batch size |

**`spatial.downstream` keys:**

| Key | Default | Description |
|-----|---------|-------------|
| `run_region_clustering` | `true` | Cluster spots by cell type composition (requires deconvolution) |
| `region_resolution` | `0.5` | Leiden resolution for region clusters |
| `region_n_neighbors` | `15` | KNN neighbours in cell-type abundance space (capped at n_cell_types − 1) |
| `run_celltype_expression` | `true` | Spearman correlation of cell type abundance vs gene expression |
| `n_marker_genes` | `20` | Top N correlated genes stored per cell type |
| `run_celltype_svg` | `true` | Moran's I on per-cell-type enriched spot subsets |
| `svg_n_genes` | `null` | Genes to test per subset (`null` = all HVGs) |
| `run_co_occurrence` | `true` | `sq.gr.co_occurrence` across distance intervals |
| `run_nhood_enrichment` | `true` | `sq.gr.nhood_enrichment` permutation test |
| `n_perms_nhood` | `1000` | Permutations for neighbourhood enrichment |
| `run_ligrec` | `true` | `sq.gr.ligrec` via OmniPath (requires internet) |
| `ligrec_n_perms` | `1000` | Permutations for ligrec |
| `ligrec_organism` | `"human"` | `"human"` or `"mouse"` |
| `run_svg_gsea` | `true` | `gseapy.prerank` on Moran's I ranked gene list |
| `svg_gsea_gene_sets` | `"GO_Biological_Process_2023"` | Enrichr gene set collection name or path to .gmt |
| `svg_gsea_organism` | `"Human"` | `"Human"` or `"Mouse"` |
| `dominant_celltype_key` | `"dominant_cell_type"` | obs column with dominant cell type per spot |
| `n_jobs` | `4` | Parallel workers for analyses that support it |

---

## 7. Tests

### `test_spatial_ingest_qc.py`
**Location:** `tests/test_spatial_ingest_qc.py`
**Covers:** `spatial_ingest` + `spatial_qc`

Key test groups:
- Auto-detection fingerprints for all 7 spatial types
- Benchmark loading via squidpy
- h5ad loading contract (layers["counts"], ENSEMBL swap, MT stripping)
- QC metric computation
- Per-threshold removal counts
- `filter_spots=True/False` behaviour
- `inplace=True/False`
- Input validation (missing obsm["spatial"], empty obs/vars)
- Summary statistics correctness

---

### `test_spatial_reduce.py`
**Location:** `tests/test_spatial_reduce.py`
**Covers:** `spatial_reduce` + `generate_spatial_reduce_report`

Key test groups:
- Contract (return type, provenance keys)
- `inplace=True/False`
- AnnData output state (all 10 expected keys present)
- Normalization logic: raw path vs benchmark path vs `normalize_total=False`
- Parameter passthrough (`n_top_genes`, `n_comps`, `flavor`)
- `n_comps` capped at min(n_obs, n_vars) - 1
- Input validation (4 error paths)
- Provenance value correctness
- Report: file generated, valid HTML, contains dataset_id and HVG count, benchmark normalization-skipped note

---

### `test_spatial_cluster.py`
**Location:** `tests/test_spatial_cluster.py`
**Covers:** `spatial_cluster` + `generate_spatial_cluster_report`

Key test groups:
- Contract (return type, provenance keys)
- `inplace=True/False`
- Clustering output state (cluster key in obs, all spots assigned, shape preserved)
- Cluster count and size consistency
- Higher resolution → more clusters (soft check)
- Annotation map: label key written, Unknown for unmapped clusters
- Custom `cluster_key`
- SVG path: `uns["moranI"]` present, DataFrame with required columns, sorted descending
- `run_svg=False`: moranI absent, no SVG provenance keys
- `svg_n_genes` caps tested genes
- Input validation (missing X_pca, missing spatial_connectivities when run_svg=True, empty obs/vars)
- Report: valid HTML, contains n_clusters, skipped-SVG section absent when run_svg=False

---

### `test_spatial_deconvolve.py`
**Location:** `tests/test_spatial_deconvolve.py`
**Covers:** `spatial_deconvolve` + `generate_spatial_deconvolve_report`
**Strategy:** cell2location training is NOT run in unit tests (too slow).
Tests cover the skip path (ref_adata=None) and a synthetic post-deconvolution state.

Key test groups:
- Skip path: all contract tests with `ref_adata=None`
- `inplace=True/False`
- Skip path: `skipped=True`, `skip_reason` present, no `q05_cell_abundance_w_sf`
- Input validation: non-AnnData, missing `layers["counts"]`, empty adata, non-AnnData ref, missing ref layer, missing cell_type_key, empty ref
- Simulated post-deconvolution state: abundance shape, non-negative values, cell type columns in obs
- `ingest._load_h5ad` contract: layers["counts"] saved, var_names swapped to ENSEMBL, MT genes stripped
- Report: valid HTML, contains dataset_id, contains cell type names, skipped notice when skipped, no provenance error

**Integration test** (skipped by default):
```python
@pytest.mark.skip(reason="requires Kuppe data + ~2h runtime")
def test_full_deconvolution_kuppe():
    ...
```
Run manually with:
```bash
python -m pytest tests/test_spatial_deconvolve.py::test_full_deconvolution_kuppe --no-header -s
```

---

### `test_spatial_downstream.py`
**Location:** `tests/test_spatial_downstream.py`
**Covers:** `spatial_downstream` + `generate_spatial_downstream_report`
**Strategy:** All analyses run on a synthetic 60-spot AnnData fixture. Slow analyses (ligrec, GSEA) are tested as individual units; co-occurrence and nhood_enrichment use squidpy on the minimal graph.

Key test groups:
- Import smoke test (module + report callable)
- `inplace=True/False` behaviour (object identity check)
- `TypeError` on non-AnnData input
- All analyses skip gracefully on bare AnnData (no deconvolution, no moranI)
- Region clustering: `obs["region_cluster"]` present, dtype is category, `n_regions ≥ 1`
- Cell-type expression: `uns["celltype_marker_genes"]` is dict, genes per cell type ≤ `n_marker_genes`
- Neighbourhood enrichment: `uns["dominant_cell_type_nhood_enrichment"]` present with `"zscore"` key
- Co-occurrence: `uns["dominant_cell_type_co_occurrence"]` present
- Cell-type SVGs: `uns["celltype_svg"]` is dict (squidpy required; skipped if absent)
- Provenance structure: `module`, `timestamp`, `params`, `analyses` all present
- Report raises `ValueError` when `uns["omicsage_spatial_downstream"]` absent
- Report generates valid HTML with expected section headings

---

## 8. AnnData State Reference

The table below shows the complete state of `adata` after each pipeline step.

| Key | After ingest | After QC | After reduce | After cluster | After deconvolve | After downstream |
|-----|:---:|:---:|:---:|:---:|:---:|:---:|
| `X` | raw counts | raw counts | normalized+log1p | normalized+log1p | normalized+log1p | normalized+log1p |
| `layers["counts"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obsm["spatial"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["spatial"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["total_counts"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["n_genes_by_counts"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["pct_counts_mt"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `obs["qc_pass"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `var["highly_variable"]` | — | — | ✅ | ✅ | ✅ | ✅ |
| `obsm["X_pca"]` | — | — | ✅ | ✅ | ✅ | ✅ |
| `obsp["spatial_connectivities"]` | — | — | ✅ | ✅ | ✅ | ✅ |
| `obs["spatial_cluster"]` | — | — | — | ✅ | ✅ | ✅ |
| `uns["moranI"]` | — | — | — | ✅ (if run_svg) | ✅ | ✅ |
| `obsm["q05_cell_abundance_w_sf"]` | — | — | — | — | ✅ (if ref) | ✅ (if ref) |
| `obs[cell_type columns]` | — | — | — | — | ✅ (if ref) | ✅ (if ref) |
| `obs["region_cluster"]` | — | — | — | — | — | ✅ (if ref) |
| `obsm["X_umap_celltype"]` | — | — | — | — | — | ✅ (if ref) |
| `uns["celltype_marker_genes"]` | — | — | — | — | — | ✅ (if ref) |
| `uns["celltype_svg"]` | — | — | — | — | — | ✅ (if ref) |
| `uns["*_nhood_enrichment"]` | — | — | — | — | — | ✅ (if ref) |
| `uns["*_co_occurrence"]` | — | — | — | — | — | ✅ (if ref) |
| `uns["*_ligrec"]` | — | — | — | — | — | ✅ (if ref) |
| `uns["svg_gsea"]` | — | — | — | — | — | ✅ |
| `uns["omicsage_spatial_ingest"]` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_qc"]` | — | ✅ | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_reduce"]` | — | — | ✅ | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_cluster"]` | — | — | — | ✅ | ✅ | ✅ |
| `uns["omicsage_spatial_deconvolve"]` | — | — | — | — | ✅ | ✅ |
| `uns["omicsage_spatial_downstream"]` | — | — | — | — | — | ✅ |

`*` = `dominant_celltype_key` prefix (default: `dominant_cell_type`).
`if ref` = requires `ref_adata` to have been passed to `spatial_deconvolve`; all downstream deconvolution-dependent analyses skip gracefully if absent.

---

## 9. Key Decisions

| Decision | Choice | Reason |
|----------|--------|--------|
| Normalization | `normalize_total` + `log1p` | SCTransform is R-only; Python/scanpy standard |
| HVG selection | `flavor="seurat"`, `n_top_genes=3000` | `"seurat"` safe for RNA (no TF-IDF overflow); 3000 standard for Visium |
| Spatial graph | `coord_type=None` (auto-detect) | Safer than hardcoding `"grid"`; squidpy auto-selects grid when `uns["spatial"]` present |
| Clustering graph | KNN from PCA, NOT spatial graph | Cluster by transcriptomic similarity; spatial graph reserved for SVGs |
| SVG p-values | `n_perms=None` (analytical) | Fast; permutation testing available but slow for 3000+ genes |
| Deconvolution | cell2location | sc-best-practices recommended; best on simulated benchmarks |
| GPU/CPU | Auto-detected by PyTorch | No explicit `use_gpu` flag; CPU fallback automatic |
| Gene IDs | ENSEMBL IDs in `var_names` after h5ad ingest | cell2location requires ENSEMBL to map between spatial and reference |
| MT genes | Stripped at ingest into `obsm["MT"]` | cell2location recommendation; MT% QC still computed from original X |
| Benchmark deconvolution | Graceful skip when `ref_adata=None` | squidpy benchmark has no paired scRNA-seq reference |
| Annotation map | Optional, no hardcoded cluster→celltype maps | Cluster numbering is non-deterministic across runs |
| `inplace` default | `False` | Safe default; avoids unexpected mutation |
| Report structure | `<header>/<main>/<footer>` with `<section>` blocks | Matches RNA pipeline exactly; required for combined report tab extraction |
| Region clustering graph | `key_added="neighbors_celltype"` in `sc.pp.neighbors` | Prevents overwriting gene-expression KNN; Leiden and UMAP for region clusters read this graph explicitly |
| UMAP key for region clusters | `obsm["X_umap_celltype"]` (not `X_umap`) | Avoids overwriting any gene-expression UMAP from earlier steps; original `X_umap` saved and restored |
| Cell-type SVG spatial graph | Recomputed per subset with `sq.gr.spatial_neighbors` | The original `obsp["spatial_connectivities"]` is invalid after spot subsetting |
| SVG GSEA p-values | `n_perms=None` (analytical) for per-subset Moran's I | Speed; running 1000 permutations × n_cell_types is prohibitive in a pipeline |
| ligrec gene symbol swap | Temporary copy with `var_names = var["feature_name"]` | OmniPath matches gene symbols; ENSEMBL IDs in `var_names` after ingest would match nothing. Copy discarded after result is transferred back |
| GSEA dtype sanitisation | `_sanitize_gsea_df()` before `uns["svg_gsea"]` assignment | gseapy returns numeric columns as `object` dtype; `write_h5ad` crashes trying to serialise floats as HDF5 variable-length strings |
| Downstream predecessor | `"cluster"` (minimum); `resolve_input` upgrades to deconvolve if available | Allows downstream to run on any dataset, even without a scRNA-seq reference; Tier 2 analyses skip gracefully |
