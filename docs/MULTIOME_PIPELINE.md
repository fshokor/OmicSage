# OmicSage — Multiome Pipeline Documentation

> Phase 5–6 — RNA + ATAC joint analysis + GRN inference  
> Last updated: June 2026  
> GitHub: https://github.com/fshokor/OmicSage

---

## Overview

The multiome pipeline processes paired single-cell RNA + ATAC data from a single experiment where both modalities are measured in the same cells (e.g. 10x Genomics Multiome). Because barcodes are shared 1:1, label transfer between modalities is exact and no probabilistic matching is needed.

The pipeline runs six sequential steps:

```
atac_qc → atac_reduce → atac_annotate → multiome_integration → multiome_deg → multiome_grn
```

Each step reads from the previous step's checkpoint, writes its own checkpoint, and generates a per-step HTML report. A combined tabbed report is assembled at the end.

**Primary benchmark dataset:** GSE194122 — NeurIPS 2021 BMMC multiome (RNA + ATAC, 13 donors, 69k cells).

---

## File Map

```
pipeline/modules/multiome/
    atac_qc.py                    Step 1 — ATAC quality control
    atac_reduce.py                Step 2 — TF-IDF → LSI → UMAP → Leiden
    atac_annotate.py              Step 3 — gene activity + label transfer
    multiome_integration.py       Step 4 — MOFA+ / MultiVI joint embedding
    multiome_deg.py               Step 5 — RNA DEG + differential accessibility
    multiome_grn.py               Step 6 — GRN inference (pyscenic + decoupler)

reports/templates/multiome/
    atac_qc_report.py             Per-step HTML report for Step 1
    atac_reduce_report.py         Per-step HTML report for Step 2
    atac_annotate_report.py       Per-step HTML report for Step 3
    multiome_integration_report.py Per-step HTML report for Step 4
    multiome_deg_report.py        Per-step HTML report for Step 5
    multiome_grn_report.py        Per-step HTML report for Step 6
    multiome_combined_report.py   Tabbed combined report (all steps)

tests/
    test_atac_qc.py
    test_atac_reduce.py
    test_atac_annotate.py
    test_multiome_integration.py
    test_multiome_deg.py
    test_multiome_grn.py
    test_run_multiome_pipeline.py

run_multiome_pipeline.py          CLI pipeline runner
config/runs/GSE194122_multiome.yaml  Dataset config
subset_multiome_raw.py            Utility — subset raw data to N batches
```

---

## Quick Start

```bash
conda activate omicsage
cd ~/OmicSage

# (Optional) subset to 6 batches first if memory is limited
python subset_multiome_raw.py

# Run the full pipeline
python run_multiome_pipeline.py --config config/runs/GSE194122_multiome.yaml

# Run a single step
python run_multiome_pipeline.py --config config/runs/GSE194122_multiome.yaml --step atac_reduce

# Run a range of steps
python run_multiome_pipeline.py --config config/runs/GSE194122_multiome.yaml \
    --from-step atac_reduce --to-step atac_annotate

# Force re-run ignoring cache
python run_multiome_pipeline.py --config config/runs/GSE194122_multiome.yaml \
    --step atac_qc --force
```

---

## Data Preparation — `subset_multiome_raw.py`

**Location:** `subset_multiome_raw.py` (repo root)

Standalone utility that subsets the raw multiome h5ad to a selected set of batches. Produces a single output file with identical structure to the input — same obs columns, same var columns, same `.X` — just fewer cells. Intended for development and testing on machines with limited RAM.

**Input:** `data/benchmark/GSE194122_multiome_raw_only.h5ad`  
**Output:** `data/benchmark/GSE194122_multiome_raw_only_6batch.h5ad`

**Usage:**

```bash
# Default — 6 largest batches
python subset_multiome_raw.py

# Custom batches
python subset_multiome_raw.py --batches s4d8 s4d1 s3d10
```

**Default batches selected** (top 6 by cell count):

| Batch | Cells |
|-------|-------|
| s4d8  | 9,876 |
| s4d1  | 8,023 |
| s3d10 | 6,781 |
| s1d2  | 6,740 |
| s1d1  | 6,224 |
| s2d4  | 6,111 |
| **Total** | **43,755** |

---

## Step 1 — `atac_qc.py`

**Location:** `pipeline/modules/multiome/atac_qc.py`  
**Input:** Raw ATAC AnnData (`.X` = raw peak counts)  
**Output:** `multiome_01_qc_atac.h5ad`  
**Report:** `multiome_01_qc_report.html`

Performs ATAC-specific quality control on the peak count matrix. Cell filtering is **disabled by default** (`filter_cells=False`) because the RNA QC step has already removed low-quality barcodes — this step computes and records ATAC metrics only.

### What it does

1. Renames CellRanger-ARC obs columns to OmicSage-namespaced equivalents
2. Preserves ground-truth labels as `obs["cell_type_groundtruth"]`
3. Computes `n_peaks_by_counts` via `sc.pp.calculate_qc_metrics`
4. Optionally runs Scrublet on the peak matrix → `obs["atac_predicted_doublet"]`
5. Applies QC threshold flags (does not drop cells unless `filter_cells=True`)
6. Filters features: keeps peaks present in ≥ `min_cells` cells
7. Saves raw counts to `layers["counts"]`
8. Writes provenance to `uns["omicsage_atac_qc"]`

### CellRanger-ARC column mapping

| Input column | OmicSage column |
|---|---|
| `ATAC_nCount_peaks` | `total_peak_counts` |
| `ATAC_atac_fragments` | `total_fragment_counts` |
| `ATAC_reads_in_peaks_frac` | `reads_in_peaks_frac` |
| `ATAC_blacklist_fraction` | `blacklist_fraction` |
| `ATAC_nucleosome_signal` | `nucleosome_signal` |

### Default thresholds (NeurIPS 2021 BMMC)

| Parameter | Default | Source |
|---|---|---|
| `min_peaks` | 750 | sc-best-practices ATAC chapter |
| `max_peaks` | 500,000 | sc-best-practices ATAC chapter |
| `min_peak_counts` | 1,500 | sc-best-practices ATAC chapter |
| `max_peak_counts` | 100,000 | sc-best-practices ATAC chapter |
| `max_nucleosome_signal` | 2.0 | sc-best-practices ATAC chapter |
| `min_cells` | 15 | feature-level filter |

### API

```python
from pipeline.modules.multiome.atac_qc import atac_qc

atac_qcd, metrics = atac_qc(
    adata,
    min_peaks=750,
    max_peaks=500_000,
    min_peak_counts=1_500,
    max_peak_counts=100_000,
    max_nucleosome_signal=2.0,
    min_cells=15,
    run_scrublet=True,
    filter_cells=False,   # True to actually drop low-quality cells
    inplace=False,
)
```

### Outputs written

| Key | Description |
|---|---|
| `obs["total_peak_counts"]` | Total peak counts per cell |
| `obs["nucleosome_signal"]` | Nucleosome signal (from CellRanger-ARC) |
| `obs["atac_predicted_doublet"]` | Scrublet doublet flag (bool) |
| `obs["cell_type_groundtruth"]` | NeurIPS ground truth label |
| `obs["pass_qc"]` | True if cell passes all thresholds |
| `layers["counts"]` | Raw peak counts (copy of input .X) |
| `uns["omicsage_atac_qc"]` | Provenance block |

> **Memory note:** Scrublet runs PCA on the full peak matrix and is the most memory-intensive part of this step. Set `run_scrublet: false` in the config if hitting memory limits.

---

## Step 2 — `atac_reduce.py`

**Location:** `pipeline/modules/multiome/atac_reduce.py`  
**Input:** `multiome_01_qc_atac.h5ad`  
**Output:** `multiome_02_reduce_atac.h5ad`  
**Report:** `multiome_02_reduce_report.html`

Computes the ATAC dimensionality reduction. LSI (Latent Semantic Indexing) is the standard method for scATAC-seq, analogous to PCA for scRNA-seq.

### What it does

1. **TF-IDF normalisation** — normalises across cells (sequencing depth) and peaks (peak accessibility) using the formula `log(TF + 1) × log(1 + N / (1 + df))`
2. **LSI via TruncatedSVD** — `sklearn.decomposition.TruncatedSVD` on the TF-IDF matrix
3. **Drop component 1** — LSI component 1 always correlates with sequencing depth, not biology; dropped unconditionally per Signac, sc-best-practices, and Seurat v5
4. **L2 row normalisation** on LSI components (standard for scATAC)
5. **Neighbor graph** on `X_lsi` → `sc.pp.neighbors`
6. **UMAP** → `obsm["X_umap_atac"]` (namespaced to avoid collision with RNA)
7. **Leiden clustering** → `obs["atac_leiden"]` (namespaced)

### API

```python
from pipeline.modules.multiome.atac_reduce import atac_reduce

atac_reduced, metrics = atac_reduce(
    adata,
    n_components=50,         # components 2–50 used (component 1 dropped)
    n_neighbors=15,
    leiden_resolution=0.5,
    use_raw_counts=True,     # use layers["counts"] as TF-IDF input
    random_state=0,
    inplace=False,
)
```

### Outputs written

| Key | Description |
|---|---|
| `layers["tf_idf"]` | TF-IDF normalised peak matrix |
| `obsm["X_lsi"]` | LSI embedding, components 2–N (n_components − 1 dims) |
| `obsm["X_umap_atac"]` | UMAP from LSI neighbors |
| `obs["atac_leiden"]` | Leiden cluster labels |
| `uns["omicsage_atac_reduce"]` | Provenance block |

### Key design decisions

**Why LSI not PCA?** ATAC peak matrices are binary/sparse; TF-IDF + SVD (LSI) is the established ATAC normalisation and dimensionality reduction approach (Cusanovich et al. 2015, Signac, sc-best-practices).

**Why drop component 1?** Verified in every major ATAC tool (Signac DepthCor, sc-best-practices, Seurat v5 `RunSVD dims=2:30`). It captures library size variation, not chromatin biology. This is unconditional and not configurable.

**Why namespaced keys?** When RNA and ATAC are combined into a MuData for integration, `X_umap_atac` and `atac_leiden` must not overwrite `X_umap` and `leiden` from the RNA pipeline.

---

## Step 3 — `atac_annotate.py`

**Location:** `pipeline/modules/multiome/atac_annotate.py`  
**Input:** `multiome_02_reduce_atac.h5ad` + `rna_input` (annotated RNA h5ad)  
**Output:** `multiome_03_annotate_atac.h5ad`  
**Report:** `multiome_03_annotate_report.html`

Annotates ATAC cells with two complementary outputs: a proxy gene expression matrix derived from chromatin accessibility, and cell type labels transferred from the RNA modality.

### Step A — Gene activity scores

Converts the peak × cell count matrix into a proxy gene expression matrix by summing peak counts that overlap each gene body + a configurable upstream promoter window.

**Peak-to-gene mapping priority:**
1. `peak_annotation=` argument (explicit override)
2. `uns["atac"]["peak_annotation"]` (populated by muon from 10x h5 files)
3. Coordinate-based proximity grouping (fallback — emits `UserWarning`; gene names will be synthetic `region_chr_start_end` labels)

> **Memory note:** The gene activity computation keeps the peak matrix sparse throughout. Individual peak columns are sliced and summed in CSC format — the full matrix is never densified.

### Step B — Label transfer

For multiome data, RNA and ATAC barcodes are 1:1. The transfer uses majority vote per ATAC Leiden cluster:

1. For each `atac_leiden` cluster, find all barcodes that also appear in the RNA AnnData
2. Collect their `obs["cell_type_vote"]` labels
3. Assign the majority label to all cells in the cluster
4. Clusters with no RNA barcode overlap → `"Unknown"`

This is simpler and more reliable than KNN label transfer (the Seurat/Signac approach for **unpaired** data) because multiome barcodes match exactly.

### API

```python
from pipeline.modules.multiome.atac_annotate import annotate_atac

atac_annotated, metrics = annotate_atac(
    atac,
    rna=rna,                          # RNA AnnData with obs["cell_type_vote"]
    peak_annotation=None,             # override uns["atac"]["peak_annotation"]
    promoter_upstream_bp=2000,
    min_peaks_per_gene=1,
    leiden_key="atac_leiden",
    rna_label_key="cell_type_vote",
    inplace=False,
)
```

### Outputs written

| Key | Description |
|---|---|
| `obsm["gene_activity"]` | Float32 array (n_cells × n_genes) — proxy gene expression |
| `uns["gene_activity_var_names"]` | Gene names for gene_activity columns |
| `obs["atac_celltype"]` | Transferred cell type label (or "Unknown") |
| `uns["omicsage_atac_annotate"]` | Provenance block |

### Namespace conventions

| obs column | Written by | Meaning |
|---|---|---|
| `atac_celltype` | `atac_annotate` | Transferred RNA label |
| `cell_type_vote` | RNA pipeline | Never written here — collision guard |
| `cell_type_groundtruth` | `atac_qc` | NeurIPS ground truth |
| `atac_leiden` | `atac_reduce` | ATAC Leiden clusters |

---

## Step 4 — `multiome_integration.py`

**Location:** `pipeline/modules/multiome/multiome_integration.py`  
**Input:** `multiome_03_annotate_atac.h5ad` + `rna_input`  
**Output:** `multiome_04_integration.h5mu`  
**Report:** `multiome_04_integration_report.html`

Assembles RNA and ATAC into a `MuData` and computes a joint low-dimensional embedding. Two methods are available.

### Method A — MOFA+ (default)

Multi-Omics Factor Analysis. Linear factor model; fast, interpretable, batch-aware via `groups_label`. Uses `mu.tl.mofa` from muon.

- **Input RNA:** `mdata["rna"].X` (log1p-normalised)
- **Input ATAC:** `mdata["atac"].X` (TF-IDF)
- **Output:** `obsm["X_mofa"]`, `obsm["X_umap_mofa"]`
- **Reference:** Argelaguet et al., Genome Biology 2020

### Method B — MultiVI

Deep generative model for paired RNA + ATAC. Models RNA with negative-binomial, ATAC with Bernoulli (peaks binarised). Requires scvi-tools.

- **Input RNA:** `mdata["rna"].layers["counts"]` (raw counts — NB model)
- **Input ATAC:** `mdata["atac"].layers["counts"]` (raw counts — binarised internally)
- **Setup:** `MULTIVI.setup_mudata` with `modalities=` dict routing each layer to the correct sub-AnnData
- **Output:** `obsm["X_multivi"]`, `obsm["X_umap_multivi"]`
- **Reference:** Ashuach et al., Nature Methods 2023

### API

```python
from pipeline.modules.multiome.multiome_integration import run_mofa, run_multivi

# MOFA+
mdata, metrics = run_mofa(
    mdata,
    batch_key="batch",
    n_factors=15,
    random_state=0,
    inplace=False,
)

# MultiVI
mdata, metrics = run_multivi(
    mdata,
    batch_key="batch",
    n_latent=20,
    max_epochs=500,
    random_state=0,
    inplace=False,
)
```

### Outputs written

| Method | Embedding key | UMAP key | Provenance key |
|---|---|---|---|
| MOFA+ | `obsm["X_mofa"]` | `obsm["X_umap_mofa"]` | `uns["omicsage_mofa"]` |
| MultiVI | `obsm["X_multivi"]` | `obsm["X_umap_multivi"]` | `uns["omicsage_multivi"]` |

> **Note:** `obsm["X_umap"]` is always renamed to the namespaced key immediately after `sc.tl.umap` — it is never left on the MuData.

### Batch key handling

`mdata.obs[batch_key]` is pushed from `mdata["rna"].obs[batch_key]` before calling MOFA+ or MultiVI. This is required because muon namespaces modality obs columns as `rna:batch` / `atac:batch` at the top level — neither muon nor scvi can see the batch column without this explicit push.

---

## Step 5 — `multiome_deg.py`

**Location:** `pipeline/modules/multiome/multiome_deg.py`  
**Input:** `multiome_04_integration.h5mu`  
**Output:** `multiome_05_deg.h5mu`  
**Report:** `multiome_05_deg_report.html`

Runs two complementary differential analyses grouped by ATAC-defined cell type labels.

### Step A — RNA DEG

Wilcoxon rank-sum test on RNA, grouped by `atac_celltype` labels transferred in Step 3. Identifies transcriptional programs underlying each chromatin-accessibility-defined population.

- Uses `mdata["rna"].layers["logcounts"]` (warns and falls back to `.X` if absent)
- Transfers `obs["atac_celltype"]` from ATAC to RNA via shared barcodes
- Results stored in `mdata["rna"].uns["rank_genes_groups_rna"]`
- Returns `dict[cell_type → DataFrame(gene, score, pval, logfc, pval_adj)]`

### Step B — DCA (Differential Chromatin Accessibility)

Wilcoxon rank-sum test on ATAC peak counts, same grouping. Identifies which genomic regions are differentially open per cell type.

- Uses `atac.layers["counts"]` raw counts (warns and falls back to TF-IDF `.X` if absent)
- Raw counts swapped into `.X` for the test, then restored in `try/finally`
- Results stored in `mdata["atac"].uns["rank_genes_groups_dca"]`
- Returns `dict[cell_type → DataFrame(peak, score, pval, logfc, pval_adj)]`
- Feature column named `peak` (not `gene`) to distinguish from RNA DEG

### Groupby priority

`atac_celltype` → `atac_leiden` → `ValueError`

### API

```python
from pipeline.modules.multiome.multiome_deg import multiome_deg

mdata, deg_dict = multiome_deg(
    mdata,
    groupby="atac_celltype",
    leiden_fallback="atac_leiden",
    method="wilcoxon",
    min_logfc=0.25,
    max_pval_adj=0.05,
    n_genes=200,
    exclude_gene_prefixes=["MT-", "RPL", "RPS"],
    exclude_peak_prefixes=["chrM"],
    inplace=False,
)
```

### Return dict keys

| Key | Type | Description |
|---|---|---|
| `rna_deg` | `dict[str, DataFrame]` | RNA DEG per cell type (`gene` column) |
| `dca` | `dict[str, DataFrame]` | DCA per cell type (`peak` column) |
| `rna_summary` | `DataFrame` | Top 5 genes per cell type (long format) |
| `dca_summary` | `DataFrame` | Top 5 peaks per cell type (long format) |
| `provenance` | `dict` | Same as `uns["omicsage_multiome_deg"]` |
| `input_type` | `str` | `"mudata"` or `"anndata"` |

### AnnData fallback

When a bare ATAC `AnnData` is passed instead of `MuData`, RNA DEG is skipped with a `UserWarning` and only DCA is computed.

---

## Step 6 — `multiome_grn.py`

**Location:** `pipeline/modules/multiome/multiome_grn.py`  
**Input:** `multiome_05_deg.h5mu`  
**Output:** `multiome_06_grn.h5mu`  
**Report:** `multiome_06_grn_report.html`

Infers gene regulatory networks (GRNs) from the paired RNA + ATAC data using a
`pyscenic` + `decoupler` stack. Full SCENIC+ is deferred (see Known Limitations).

### What it does

1. **Fetches CollecTRI TF–target network** from the OmniPath REST API at runtime
   (~47k edges, single-gene TFs only — complex entries like `MYC_MAX` excluded)
2. **B1 — ATAC motif enrichment:** extracts top `n_top_peaks` DCA peaks per cell
   type from `atac.uns["rank_genes_groups_dca"]`; scores cells for TF activity
   via `decoupler.mt.aucell` on the accessible peak set
3. **B2 — RNA TF activity:** builds TF regulons via pyscenic correlation fallback
   (GRNBoost2 if `arboreto` is installed); scores cells via `pyscenic.aucell`
4. **B3 — GRN table:** merges RNA and ATAC scores into a unified edge DataFrame
   (TF → target gene, per cell type) sorted by combined score

### API

```python
from pipeline.modules.multiome.multiome_grn import multiome_grn

mdata, grn_dict = multiome_grn(
    mdata,
    deg_dict=deg_dict,          # provenance dict from multiome_deg
    motif_db="jaspar",          # "jaspar" = OmniPath CollecTRI; or path to .feather
    groupby="atac_celltype",
    n_top_peaks=500,
    min_cells=10,
    random_state=0,
    inplace=False,
)
```

### Outputs written

| Key | Description |
|---|---|
| `obsm["X_aucell_rna"]` | pyscenic AUCell scores (n_cells × n_RNA_regulons) |
| `obsm["X_aucell_atac"]` | decoupler AUCell scores (n_cells × n_ATAC_TFs) |
| `uns["grn_network"]` | GRN edge table as dict-of-lists (JSON-serialisable) |
| `uns["omicsage_grn"]` | Provenance block |

### Return dict keys

| Key | Type | Description |
|---|---|---|
| `provenance` | `dict` | Same as `uns["omicsage_grn"]` |
| `n_tfs_rna` | `int` | Number of RNA regulons scored |
| `n_tfs_atac` | `int` | Number of ATAC motif TFs scored |
| `n_grn_edges` | `int` | Total rows in GRN edge table |
| `grn_df` | `DataFrame` | Columns: tf, target_gene, rna_score, atac_score, combined_score, cell_type |

### Dependencies

```bash
pip install "pyscenic>=0.12.1" "ctxcore>=0.2.0" "decoupler>=2.0.0"
```

All three are optional at runtime — if any are missing the function returns a
valid empty result with a `UserWarning` rather than raising.

### CollecTRI fetch

CollecTRI is downloaded at runtime from OmniPath with no local database required:

```
https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1
```

Requires internet access. The network is fetched once per `multiome_grn` call
and is not cached to disk. Offline use: pass a pre-downloaded `.feather` file
path as `motif_db`.

### Known limitations

| Limitation | Detail |
|---|---|
| ATAC scores are flat (all 1.0) | No coordinate-level motif–peak overlap; all DCA peaks treated as equally bound by all TFs. Requires cisTarget DB or `pybedtools` intersection to resolve. |
| RNA TFs = 0 on some datasets | pyscenic AUCell requires TFs present in `var_names`; depends on CollecTRI overlap with HVG subset (2000 genes). Larger gene sets improve recall. |
| No cisTarget database | Full SCENIC+ motif enrichment requires a 3–10 GB feather file not included. Pass path as `motif_db` param when available. |
| GRNBoost2 not used by default | `arboreto` not in base omicsage env. Install separately for higher-quality RNA regulons: `pip install arboreto` |

### Upgrade path to full SCENIC+

```bash
# Separate env recommended to avoid dependency conflicts
conda create -n omicsage_grn python=3.11
pip install pycistarget scenicplus
# Download cisTarget JASPAR feather (~3 GB) from:
# https://resources.aertslab.org/cistarget/databases/
# Then pass to multiome_grn:
mdata, grn = multiome_grn(mdata, deg_dict, motif_db="path/to/cistarget.feather")
```

---

## Pipeline Runner — `run_multiome_pipeline.py`

**Location:** `run_multiome_pipeline.py` (repo root)

### Step registry

| Step | Checkpoint | Format |
|---|---|---|
| `atac_qc` | `multiome_01_qc_atac.h5ad` | AnnData |
| `atac_reduce` | `multiome_02_reduce_atac.h5ad` | AnnData |
| `atac_annotate` | `multiome_03_annotate_atac.h5ad` | AnnData |
| `multiome_integration` | `multiome_04_integration.h5mu` | MuData |
| `multiome_deg` | `multiome_05_deg.h5mu` | MuData |
| `multiome_grn` | `multiome_06_grn.h5mu` | MuData |

### CLI flags

| Flag | Description |
|---|---|
| `--config` | Path to YAML config (required) |
| `--step` | Run exactly one step |
| `--from-step` | Start from this step (inclusive) |
| `--to-step` | Stop at this step (inclusive) |
| `--force` | Re-run even if cached output exists |

### Config structure

```yaml
dataset:
  id:   GSE194122_multiome
  name: "BMMC Multiome (NeurIPS 2021)"

paths:
  atac_input:    data/processed/GSE194122/01_qc_atac.h5ad
  rna_input:     data/processed/GSE194122/05_annotated.h5ad
  processed_dir: data/processed/GSE194122_multiome
  reports_dir:   reports/GSE194122_multiome

multiome:
  batch_key: batch

steps:
  atac_qc:
    enabled: true
    params:
      min_peaks: 750
      run_scrublet: false    # disable to save memory on large datasets

  atac_reduce:
    enabled: true
    params:
      n_components: 50
      leiden_resolution: 0.5

  atac_annotate:
    enabled: true
    params:
      promoter_upstream_bp: 2000
      rna_label_key: cell_type_vote

  multiome_integration:
    enabled: true
    params:
      method: mofa          # mofa | multivi
      batch_key: batch
      n_factors: 15

  multiome_grn:
    enabled: true
    params:
      motif_db: jaspar          # "jaspar" = OmniPath CollecTRI; or path to .feather
      groupby: atac_celltype
      n_top_peaks: 500
      min_cells: 10
```

### Two special steps

`atac_annotate` and `multiome_integration` both require `paths.rna_input` in addition to their predecessor checkpoint. The runner validates this before execution and raises a clear error if the file is missing.

---

## Combined Report — `multiome_combined_report.py`

**Location:** `reports/templates/multiome/multiome_combined_report.py`  
**Output:** `multiome_00_combined_report.html`

Assembles all available per-step HTML reports into a single self-contained tabbed file. Called automatically at the end of `run_multiome_pipeline.py`. Can also be run standalone to rebuild from existing reports.

### Tab order

| Tab | Step report file |
|---|---|
| 🔬 ATAC QC | `multiome_01_qc_report.html` |
| 🔭 Reduce ATAC | `multiome_02_reduce_report.html` |
| 🏷️ Annotate ATAC | `multiome_03_annotate_report.html` |
| 🔗 Integration | `multiome_04_integration_report.html` |
| 📊 DEG / DCA | `multiome_05_deg_report.html` |
| 🧬 GRN | `multiome_06_grn_report.html` |

Only tabs for reports that exist on disk are rendered — partial pipeline runs produce a partial combined report without errors.

```bash
# Rebuild standalone
python -m reports.templates.multiome.multiome_combined_report \
    --reports-dir reports/GSE194122_multiome \
    --dataset-name "BMMC Multiome (NeurIPS 2021)"
```

---

## Namespace Conventions

All multiome keys are namespaced to avoid collision with the RNA pipeline when both are combined into a MuData.

| Key | Written by | Collides with |
|---|---|---|
| `obsm["X_lsi"]` | `atac_reduce` | (no RNA equivalent) |
| `obsm["X_umap_atac"]` | `atac_reduce` | RNA `obsm["X_umap"]` |
| `obs["atac_leiden"]` | `atac_reduce` | RNA `obs["leiden"]` |
| `obs["atac_celltype"]` | `atac_annotate` | RNA `obs["cell_type_vote"]` |
| `obsm["gene_activity"]` | `atac_annotate` | (no RNA equivalent) |
| `obsm["X_mofa"]` | `multiome_integration` | RNA `obsm["X_pca_harmony"]` |
| `obsm["X_umap_mofa"]` | `multiome_integration` | RNA `obsm["X_umap"]` |
| `obsm["X_multivi"]` | `multiome_integration` | RNA `obsm["X_pca_harmony"]` |
| `obsm["X_umap_multivi"]` | `multiome_integration` | RNA `obsm["X_umap"]` |
| `obsm["X_aucell_rna"]` | `multiome_grn` | (no RNA pipeline equivalent) |
| `obsm["X_aucell_atac"]` | `multiome_grn` | (no RNA pipeline equivalent) |
| `uns["grn_network"]` | `multiome_grn` | (no RNA pipeline equivalent) |

---

## Checkpoint Flow

```
paths.atac_input  ─────────────────────────────────────────────────────────────────┐
                                                                                    │
                  ┌──────────────┐   h5ad   ┌──────────────┐   h5ad               │
                  │  atac_qc     │ ───────▶  │  atac_reduce │ ───────▶  ...        │
                  └──────────────┘           └──────────────┘                      │
                                                                                    │
paths.rna_input  ──────────────────────────────────────────────────────────────────┤
                                                                          (Step 3)  │
                  ┌────────────────┐   h5ad   ┌──────────────────────┐   h5mu     │
  ... ──────────▶ │  atac_annotate │ ───────▶  │ multiome_integration │ ───────▶  │
                  └────────────────┘           └──────────────────────┘           │
                  ↑ also needs rna_input        ↑ also needs rna_input             │
                                                                                    │
                  ┌───────────────┐   h5mu   ┌───────────────┐   h5mu             │
  ... ──────────▶ │ multiome_deg  │ ───────▶  │ multiome_grn  │ ───────▶  done    │
                  └───────────────┘           └───────────────┘                    │
```

---

## Known Issues

| Issue | Workaround |
|---|---|
| `_sum_peaks_to_genes` OOM on large datasets | Fixed — uses CSC sparse slicing, never densifies full matrix |
| Scrublet OOM on 116k peaks | Set `run_scrublet: false` in config |
| WSL2 disconnect under memory pressure | Set `memory=12GB` in `~/.wslconfig`, restart with `wsl --shutdown` |
| WNN via `muon.pp.neighbors` hangs on small fixtures | Marked `@pytest.mark.slow`, deferred to full-data run |
| MultiVI training too slow for CI | MultiVI tests commented out — uncomment for local GPU runs |
| `peak_annotation` absent after h5ad round-trip | Load from original 10x h5 file; coordinate fallback emits warning |
| GRN ATAC scores all 1.0 | No coordinate-level motif–peak overlap without cisTarget DB; see Step 6 Known Limitations |
| GRN RNA TFs = 0 | CollecTRI TF overlap with 2000-gene HVG subset may be low; use larger gene set |
| `dc.run_aucell` AttributeError | decoupler 2.x moved AUCell to `dc.mt.aucell`; fixed in current implementation |
| OmniPath `fields` parameter rejected | Use bare `?datasets=collectri&genesymbols=1` URL; fixed in current implementation |

---

## References

- sc-best-practices ATAC QC: https://www.sc-best-practices.org/chromatin_accessibility/quality_control.html
- sc-best-practices paired integration (MultiVI): https://www.sc-best-practices.org/multimodal_integration/paired_integration.html
- sc-best-practices GRN chapter: https://www.sc-best-practices.org/chromatin_accessibility/gene_regulatory_networks_atac.html
- Signac PBMC vignette: https://stuartlab.org/signac/articles/pbmc_vignette.html
- Seurat v5 ATAC integration: https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
- MOFA+ paper: Argelaguet et al., Genome Biology 2020
- MultiVI paper: Ashuach et al., Nature Methods 2023
- SCENIC+ paper: Bravo González-Blas et al., Nature Methods 2023 — https://www.nature.com/articles/s41592-023-01938-4
- pySCENIC GitHub: https://github.com/aertslab/pySCENIC
- CollecTRI / OmniPath: https://omnipathdb.org/
- NeurIPS 2021 BMMC dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122
