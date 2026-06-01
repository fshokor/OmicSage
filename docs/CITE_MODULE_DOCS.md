# OmicSage — CITE-seq Module Reference

> Locations: `pipeline/modules/cite/` · `reports/templates/cite/` · repo root
> Phase: 4 — CITE-seq Module
> Last updated: 2026-05-26 (session 7)

This document describes every script in the CITE-seq pipeline — what it does,
what goes in, what comes out, and how it connects to the next step.

The CITE-seq pipeline processes **ADT (antibody-derived tag)** data from paired
RNA + protein experiments. It runs after the RNA pipeline QC step and operates
on the `01_qc_adt.h5ad` file that `qc.py` produces automatically when the input
dataset contains both RNA and ADT features.

**Step order**

```
normalize_adt → doublets → reduce_adt → harmony_adt → annotate_adt → integration
```

**Checkpoint files produced** (all under `data/processed/<dataset_id>/`)

| Step | Output file |
|------|-------------|
| normalize_adt | `cite_01_normalized_adt.h5ad` |
| doublets | `cite_02_doublets_adt.h5ad` |
| reduce_adt | `cite_03_reduced_adt.h5ad` |
| harmony_adt | `cite_04_harmony_adt.h5ad` |
| annotate_adt | `cite_05_annotated_adt.h5ad` |
| integration | `cite_06_integration.h5mu` (MuData) |

**Reports produced** (all under `reports/<dataset_id>/`)

| Step | Report file |
|------|-------------|
| normalize_adt | `cite_01_normalize_report.html` |
| doublets | `cite_02_doublets_report.html` |
| reduce_adt | `cite_03_reduce_report.html` |
| harmony_adt | `cite_04_harmony_report.html` |
| annotate_adt | `cite_05_annotate_report.html` |
| integration | `cite_06_integration_report.html` |

---

## 1. `adt_normalize.py`

**Location**: `pipeline/modules/cite/adt_normalize.py`

**What it does**

Normalizes raw ADT count data using **CLR (Centered Log-Ratio)** normalization.
Takes the `mdata["adt"]` AnnData from `qc.py` (raw integer counts in `.X`) and
returns a normalized AnnData with CLR values in both `.X` and `layers["adt_clr"]`,
and raw counts preserved in `layers["counts"]`.

**Why it exists**

ADT protein counts cannot be normalized the same way as RNA. Total-count
normalization (CP10K) is inappropriate because the set of proteins measured is
small and deliberately targeted — a spike in one protein genuinely changes the
relative abundance of others. CLR normalization divides each protein count by
the geometric mean of all proteins in the same cell or across cells, then
applies log1p. This is the standard approach from Seurat and sc-best-practices.

**Why CLR over DSB**

DSB (Denoised and Scaled by Background) is technically superior but requires
empty-droplet ADT counts from the unfiltered droplet matrix. Most public datasets
do not provide this. CLR is robust and well-validated on published CITE-seq
benchmarks including the NeurIPS 2021 BMMC dataset.

**Layer layout after normalize_adt()**

| Slot | Contents |
|------|----------|
| `.X` | CLR-normalized values |
| `layers["counts"]` | Raw integer counts (preserved from input) |
| `layers["adt_clr"]` | CLR-normalized values (named copy — always accessible by name) |
| `uns["omicsage_adt_normalize"]` | Provenance record with timestamp and parameters |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Raw ADT counts in `.X` — pass `mdata["adt"]` from `run_qc()` |
| `clr_axis` | int | `0` | `0` = per-protein across cells (muon default, Seurat `margin=1`); `1` = per-cell across proteins (Seurat `margin=2`) |
| `dsb_empty_adata` | AnnData | None | Empty-droplet AnnData for DSB normalization — not yet implemented |
| `isotype_controls` | list[str] | None | Isotype control antibody names for DSB — not yet implemented |
| `inplace` | bool | False | Modify input AnnData in place; default makes a copy |

**CLR axis choice**

- `axis=0` (default) — recommended by muon and sc-best-practices for bulk CITE-seq normalization.
  Pools all cells to estimate per-protein background. Use this unless you have
  a specific reason to prefer `axis=1`.
- `axis=1` — normalizes within each cell independently. Better when cells have
  very different total ADT counts (e.g. when mixing very different cell types).
  Equivalent to Seurat `NormalizeData(margin=2)`.

**Input**

```
adata  : AnnData  →  raw integer ADT counts in adata.X
                     typically mdata["adt"] from run_qc()
```

**Output**

```
adata_norm  : AnnData  →  CLR-normalized ADT
                           .X                     = CLR values
                           layers["counts"]       = raw counts
                           layers["adt_clr"]      = CLR values (named layer)
                           uns["omicsage_adt_normalize"]  = provenance
metrics     : dict     →  n_cells, n_proteins, clr_axis, clr_max, clr_min,
                           raw_median_total_counts_per_cell
```

**Usage**

```python
from pipeline.modules.cite.adt_normalize import normalize_adt

# Default — CLR per-protein across cells (recommended)
adt_norm, metrics = normalize_adt(mdata["adt"])

# Per-cell CLR (Seurat margin=2)
adt_norm, metrics = normalize_adt(mdata["adt"], clr_axis=1)

# In-place (no copy)
adt_norm, metrics = normalize_adt(mdata["adt"], inplace=True)

print(metrics["n_proteins"])    # 134
print(metrics["clr_max"])       # ~8.0
```

**Connects to**: `cite_normalize_report.py` for the report, then `adt_doublets.py`

---

## 2. `adt_doublets.py`

**Location**: `pipeline/modules/cite/adt_doublets.py`

**What it does**

Detects **heterotypic ADT doublets** — cells that appear to co-express mutually
exclusive lineage surface markers. A real single cell cannot simultaneously
express both CD3 (T-cell marker) and CD19 (B-cell marker) above background.
Co-expression above a CLR threshold is scored as evidence of a doublet.

**Why it exists**

Scrublet (used in `qc.py`) detects RNA-level doublets. Some multiplets pass
RNA-based filtering but are detectable in protein space because ADT has less
ambient noise. The ADT doublet score is a complementary, orthogonal quality
metric that catches multiplets that RNA QC misses.

This module only flags doublets by default — it does not remove them. This is
intentional: the analyst should inspect `obs["adt_doublet_score"]` in the
normalization report and decide whether to filter based on the observed distribution.
Set `filter_doublets=True` in the config to remove flagged cells.

**How the score is computed**

For each marker pair (e.g. CD3 / CD19):
1. Resolve each marker name to a column in `var_names` by prefix match (case-insensitive).
   `"CD3"` matches `"CD3-TotalSeqB"`, `"CD3_1"`, etc.
2. Flag the cell if **both** markers in the pair exceed the CLR threshold.
3. The score is the fraction of evaluated pairs where both markers are
   simultaneously expressed: `score = n_flagged_pairs / n_evaluated_pairs` (0.0–1.0).
4. A cell is called a doublet if `score > 0` (at least one pair flagged).

**Default marker pairs**

| Pair | Lineage conflict |
|------|----------------|
| CD3 / CD19 | T cell vs B cell |
| CD3 / CD14 | T cell vs Monocyte |

Both pairs are resolved by prefix match — if neither marker from a pair is
found in `var_names`, that pair is skipped and recorded in `metrics["pairs_skipped"]`.

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mdata` | MuData | required | Must contain `mdata["adt"]` with `layers["adt_clr"]` |
| `marker_pairs` | list[(str, str)] | `[("CD3","CD19"),("CD3","CD14")]` | Mutually exclusive marker pairs |
| `threshold` | float | `2.5` | CLR value above which a marker is "expressed"; sc-best-practices recommendation |
| `filter_doublets` | bool | False | Remove flagged cells from both `adt` and `rna` modalities |
| `inplace` | bool | False | Modify `mdata` in place; default makes a copy |

**Input**

```
mdata  : MuData  →  must contain mdata["adt"].layers["adt_clr"]
                    (output of normalize_adt)
```

**Output**

```
mdata   : MuData  →  updated object
                      mdata["adt"].obs["adt_doublet_score"]     = float 0–1
                      mdata["adt"].obs["adt_predicted_doublet"] = bool
                      mdata["adt"].uns["omicsage_adt_doublets"] = provenance
metrics : dict    →  n_cells_before, n_doublets_detected, pct_doublets,
                      pairs_evaluated, pairs_skipped, n_cells_after,
                      threshold, filter_doublets
```

**Usage**

```python
from pipeline.modules.cite.adt_doublets import detect_adt_doublets

# Default — flag only, do not remove
mdata, metrics = detect_adt_doublets(mdata)
print(f"{metrics['n_doublets_detected']} doublets ({metrics['pct_doublets']:.1f}%)")

# Flag and remove
mdata, metrics = detect_adt_doublets(mdata, filter_doublets=True)

# Custom marker pairs and threshold
mdata, metrics = detect_adt_doublets(
    mdata,
    marker_pairs=[("CD3", "CD19"), ("CD3", "CD56"), ("CD14", "CD19")],
    threshold=3.0,
    filter_doublets=False,
)
```

**Connects to**: `cite_doublets_report.py` for the report, then `adt_reduce.py`

---

## 3. `adt_reduce.py`

**Location**: `pipeline/modules/cite/adt_reduce.py`

**What it does**

Runs **PCA → neighbor graph → UMAP** on CLR-normalized ADT data. This is the
ADT-space equivalent of the RNA `reduce.py` step, but with protein-specific
defaults and strict naming conventions to prevent key collisions with RNA embeddings.

**Why it exists**

ADT data has far fewer features than RNA (typically 30–200 proteins vs 2000 HVGs).
PCA on ADT captures protein co-expression patterns rather than transcriptional
states. The resulting ADT PCA embedding is the input to Harmony batch correction
and Leiden clustering in the next two steps.

Without separate ADT-specific obsm keys (`X_pca_adt`, `X_umap_adt`), running
this step after the RNA pipeline would silently overwrite `obsm["X_pca"]` and
`obsm["X_umap"]`, corrupting the RNA embedding. The `_adt` suffix is enforced
unconditionally.

**Embedding keys written**

| Slot | Contents |
|------|----------|
| `obsm["X_pca_adt"]` | ADT PCA embedding (n_cells × n_comps_actual) |
| `obsm["X_umap_adt"]` | ADT UMAP embedding, pre-Harmony (n_cells × 2) |
| `uns["pca_adt"]` | PCA variance statistics |
| `uns["omicsage_adt_reduce"]` | Provenance record |

Note: `obsm["X_umap_adt"]` is overwritten by `adt_harmony.py` with the post-Harmony
UMAP. The pre-Harmony UMAP is visible in `cite_03_reduce_report.html`.

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mdata` | MuData | required | Must contain `mdata["adt"].layers["adt_clr"]` |
| `n_comps` | int | `50` | PCA components to compute — capped at `min(n_cells-1, n_vars-1)` |
| `n_pcs` | int | `20` | PCs used for the neighbor graph; sc-best-practices ch.36 default |
| `n_neighbors` | int | `15` | kNN graph k |
| `svd_solver` | str | `"arpack"` | SVD solver for PCA |
| `random_state` | int | `0` | Reproducibility seed |
| `inplace` | bool | False | Modify `mdata["adt"]` in place; default makes a copy |

**Why n_pcs=20 for ADT**

ADT panels have 30–200 proteins. Most biological variation is captured in the
first 15–25 PCs. Using all 50 PCs for the neighbor graph adds noise from the
low-variance tail. 20 PCs is the sc-best-practices recommendation (ch.36) and
is appropriate for the 134-protein GSE194122 panel.

**Input**

```
mdata  : MuData  →  mdata["adt"].layers["adt_clr"]  (output of normalize_adt)
```

**Output**

```
adata_reduced  : AnnData  →  mdata["adt"] with PCA and UMAP computed
                              obsm["X_pca_adt"]   = ADT PCA
                              obsm["X_umap_adt"]  = ADT UMAP (pre-Harmony)
metrics        : dict     →  n_cells, n_vars, n_comps_actual, n_pcs_used,
                              variance_explained_total, umap_computed
```

**Usage**

```python
from pipeline.modules.cite.adt_reduce import reduce_adt

adt_reduced, metrics = reduce_adt(mdata)

# Custom parameters
adt_reduced, metrics = reduce_adt(
    mdata,
    n_comps=50,
    n_pcs=20,
    n_neighbors=15,
    random_state=0,
    inplace=False,
)

print(metrics["n_comps_actual"])            # actual PCs computed (may be < 50 for small data)
print(metrics["variance_explained_total"])  # total variance in computed PCs
```

**Connects to**: `cite_reduce_report.py` for the report, then `adt_harmony.py`

---

## 4. `adt_harmony.py`

**Location**: `pipeline/modules/cite/adt_harmony.py`

**What it does**

Runs **Harmony batch correction** on the ADT PCA embedding (`obsm["X_pca_adt"]`),
writes the corrected embedding to `obsm["X_pca_harmony_adt"]`, recomputes the
neighbor graph on the corrected embedding, and overwrites `obsm["X_umap_adt"]`
with a post-Harmony UMAP.

**Why it exists**

The GSE194122 BMMC dataset was collected across multiple donors and sites in the
NeurIPS 2021 challenge. Without batch correction, donor effects dominate the ADT
embedding and cells cluster by donor rather than by cell type. Harmony aligns
embeddings across batches while preserving biological variation.

The module calls harmonypy directly rather than via `sc.external.pp.harmony_integrate`
to avoid a shape bug introduced in harmonypy ≥ 2.0.0. The scanpy wrapper does
`harmony_out.Z_corr.T` which is correct for ≤ 0.0.9 (where `Z_corr` is
PCs × cells) but wrong for ≥ 2.0.0 (where it is cells × PCs). This module
detects the transpose at runtime.

**Key naming rule**

The Harmony-corrected embedding is stored as `obsm["X_pca_harmony_adt"]`, not
`obsm["X_pca_harmony"]`. The `_adt` suffix prevents collision with the RNA
Harmony embedding (`obsm["X_pca_harmony"]`) produced by the RNA pipeline's
Harmony step. Both can coexist in the same MuData without ambiguity.

**Steps performed (in order)**

1. Validate inputs — checks `obsm["X_pca_adt"]` and `obs[batch_key]` exist
2. Run harmonypy directly on `obsm["X_pca_adt"]`
3. Detect Z_corr shape and transpose if needed (harmonypy version compatibility)
4. Store corrected embedding in `obsm["X_pca_harmony_adt"]`
5. Recompute neighbor graph using `X_pca_harmony_adt`
6. Recompute UMAP → overwrite `obsm["X_umap_adt"]` with harmony UMAP

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mdata` | MuData | required | Must contain `mdata["adt"].obsm["X_pca_adt"]` and `obs[batch_key]` |
| `batch_key` | str | required | `obs` column used as the batch variable — e.g. `"batch"` or `"donor"` |
| `n_pcs` | int | `20` | Harmony-corrected PCs used to build post-correction neighbor graph |
| `n_neighbors` | int | `15` | kNN graph k |
| `random_state` | int | `0` | Harmony random seed |
| `inplace` | bool | False | Modify `mdata["adt"]` in place; default makes a copy |

**Input**

```
mdata  : MuData  →  mdata["adt"].obsm["X_pca_adt"]  (output of reduce_adt)
                    mdata["adt"].obs[batch_key]       (batch column, e.g. "batch")
```

**Output**

```
adata_harmonised  : AnnData  →  mdata["adt"] with Harmony correction applied
                                 obsm["X_pca_harmony_adt"] = Harmony-corrected PCA
                                 obsm["X_umap_adt"]        = UMAP from harmony embedding
                                 obsp["connectivities"]    = post-harmony kNN graph
                                 obsp["distances"]         = post-harmony kNN distances
                                 uns["omicsage_adt_harmony"] = provenance
metrics           : dict     →  n_cells, n_pcs_used, n_batches, batch_values,
                                 harmony_key, umap_key, random_state
```

**Usage**

```python
from pipeline.modules.cite.adt_harmony import run_harmony_adt

adt_harmonised, metrics = run_harmony_adt(mdata, batch_key="batch")

print(metrics["n_batches"])     # number of unique batches
print(metrics["batch_values"])  # ["D1", "D2", ..., "D10"]
```

**Connects to**: `cite_harmony_report.py` for the report, then `adt_annotate.py`

---

## 5. `adt_annotate.py`

**Location**: `pipeline/modules/cite/adt_annotate.py`

**What it does**

Runs **Leiden clustering** at low resolution on the post-Harmony neighbor graph,
computes `rank_genes_groups` and a dendrogram for dotplot support, and optionally
maps cluster IDs to cell type labels when an `annotation_map` is provided.

**Why it exists**

ADT-based cell type annotation is more direct than RNA annotation because surface
protein expression patterns are less noisy and more lineage-specific. CD3+/CD4+
cells are CD4 T cells. CD3+/CD8+ cells are CD8 T cells. CD19+ cells are B cells.
A low-resolution Leiden clustering (0.1) on the ADT neighbor graph typically
produces 5–10 broad immune clusters that map cleanly to known lineages.

The `annotation_map` workflow is deliberately two-step: run with
`annotation_map=null` first, inspect `cite_05_annotate_report.html`, identify
clusters by their protein expression patterns, then re-run with the map filled in.
This avoids the need for automated marker-to-label matching that can fail on
non-standard panels.

**Naming convention**

Cell type labels are stored in `obs["adt_celltype"]`, not `obs["celltype"]` or
`obs["cell_type"]`. This prevents collision with the RNA pipeline's consensus
annotation column `obs["cell_type_vote"]` when both AnnData objects are joined
into a MuData for integration.

**Steps performed (in order)**

1. Validate inputs — checks `obsp["connectivities"]` is present
2. Leiden clustering (`flavor="igraph"`, `directed=False`) on post-Harmony graph
3. `rank_genes_groups` — marker proteins per cluster (for dotplot)
4. `dendrogram` on Leiden clusters (for ordered dotplot)
5. Optional: apply `annotation_map` → write `obs["adt_celltype"]`

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mdata` | MuData | required | Must contain `mdata["adt"].obsp["connectivities"]` |
| `annotation_map` | dict[str, str] | None | Maps Leiden cluster ID strings to cell type labels, e.g. `{"0": "CD4 T", "1": "B cell"}` |
| `resolution` | float | `0.1` | Leiden resolution — lower = fewer, broader clusters |
| `n_iterations` | int | `2` | Leiden iterations (igraph flavor) |
| `random_state` | int | `0` | Leiden seed |
| `inplace` | bool | False | Modify `mdata["adt"]` in place; default makes a copy |

**Input**

```
mdata  : MuData  →  mdata["adt"].obsp["connectivities"]  (output of harmony_adt)
```

**Output**

```
adata_annotated  : AnnData  →  mdata["adt"] with clustering applied
                                obs["leiden"]                 = Leiden cluster IDs (always)
                                obs["adt_celltype"]           = cell type labels (if annotation_map)
                                uns["rank_genes_groups"]      = marker protein results
                                uns["dendrogram_leiden"]      = cluster dendrogram
                                uns["omicsage_adt_annotate"]  = provenance
metrics          : dict     →  n_cells, n_clusters, cluster_sizes, annotated,
                                leiden_key, celltype_key, resolution
```

**Usage**

```python
from pipeline.modules.cite.adt_annotate import annotate_adt

# First run — no annotation map, just clustering
adt_annotated, metrics = annotate_adt(mdata)
# → inspect cite_05_annotate_report.html
# → identify clusters from the UMAP and cluster sizes table

# Second run — apply annotation map
annotation_map = {
    "0": "CD4 T",
    "1": "CD8 T",
    "2": "B cell",
    "3": "NK",
    "4": "Monocyte",
    "5": "DC",
}
adt_annotated, metrics = annotate_adt(mdata, annotation_map=annotation_map)
print(metrics["n_clusters"])  # 6
```

**Connecting annotation_map to the config**

After identifying clusters, add the map to `config/runs/GSE194122_cite.yaml`
under `steps.annotate_adt.params.annotation_map` and re-run that step only:

```bash
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --step annotate_adt --force
```

**Connects to**: `cite_annotate_report.py` for the report, then `cite_integration.py`

---

## 6. `cite_integration.py`

**Location**: `pipeline/modules/cite/cite_integration.py`

**What it does**

Implements two multi-modal integration methods for paired RNA + ADT data:

- **MOFA+** (Multi-Omics Factor Analysis) — a linear factor model that
  decomposes the joint RNA + ADT variance into shared latent factors.
  Fast, interpretable, no GPU required.

- **totalVI** (Total Variational Inference) — a deep generative model that
  explicitly models RNA with a negative-binomial distribution and ADT as a
  foreground/background NB mixture. Higher benchmark scores but slower and
  requires more RAM.

Both methods produce a joint latent embedding stored in `mdata.obsm`, plus a
UMAP computed from that embedding for visualization.

**Why it exists**

RNA and ADT measure complementary aspects of cell state. RNA captures
transcriptional programmes; ADT captures surface protein abundance. Analyzing
them independently misses the joint biological signal. Integration produces
a single embedding that is informed by both modalities, enabling more accurate
cell type identification than either modality alone — particularly for rare
cell types that are well-separated by protein but noisy in RNA.

**Why MOFA+ is the default**

MOFA+ runs in minutes on the full GSE194122 dataset without a GPU. The latent
factors are interpretable (you can ask "which proteins load onto factor 3?").
It handles batch effects via `groups_label=batch_key`. For most exploratory
analyses, MOFA+ produces results that are good enough to proceed.

totalVI should be preferred for a final analysis or publication because it
achieves higher scIB benchmark scores and explicitly models ADT background noise.

**Important implementation detail — MOFA+ batch key**

muon namespaces conflicting obs columns when building `mdata.obs`. If both
`mdata["rna"].obs` and `mdata["adt"].obs` have a `"batch"` column, muon
renames them `"rna:batch"` and `"adt:batch"` in `mdata.obs`. MOFA+ reads
`groups_label` from `mdata.obs` directly and cannot find `"batch"`. This
module pushes `batch_key` to `mdata.obs` before calling `mu.tl.mofa` as a
workaround. The column is left there intentionally as it is useful for
downstream visualization.

**Methods**

### run_mofa

```python
run_mofa(mdata, batch_key, n_factors=15, use_layer=None,
         random_state=0, inplace=False)
→ (MuData, dict)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mdata` | MuData | required | RNA in `mdata["rna"].X` (log1p); ADT in `mdata["adt"].X` (CLR) |
| `batch_key` | str | required | `obs` column for MOFA+ `groups_label` — typically `"batch"` |
| `n_factors` | int | `15` | Number of latent factors |
| `use_layer` | str | None | Use this layer from each modality instead of `.X` |
| `random_state` | int | `0` | MOFA+ seed |
| `inplace` | bool | False | Modify `mdata` in place |

Output keys written to `mdata`:
```
obsm["X_mofa"]        →  MOFA+ latent factors (n_cells × n_factors)
obsm["X_umap_mofa"]   →  UMAP from MOFA+ embedding (n_cells × 2)
uns["omicsage_mofa"]  →  provenance
```

### run_totalvi

```python
run_totalvi(mdata, batch_key, max_epochs=400, random_state=0, inplace=False)
→ (MuData, dict)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mdata` | MuData | required | RNA raw counts in `mdata["rna"].layers["counts"]`; ADT raw counts in `mdata["adt"].layers["counts"]`; batch column on RNA modality |
| `batch_key` | str | required | `obs` column for batch — typically `"batch"` |
| `max_epochs` | int | `400` | Training epochs — use `10` for a quick test |
| `random_state` | int | `0` | scvi-tools seed |
| `inplace` | bool | False | Modify `mdata` in place |

Output keys written to `mdata`:
```
obsm["X_totalVI"]       →  totalVI latent representation (n_cells × n_latent)
obsm["X_umap_totalVI"]  →  UMAP from totalVI embedding (n_cells × 2)
uns["omicsage_totalVI"] →  provenance
```

**Note on WNN (deferred)**

WNN (Weighted Nearest Neighbor) via `muon.pp.neighbors` is not implemented
in the current release. It hangs on small fixtures due to pynndescent behavior.
When implemented, outputs will be:
```
mdata.obsm["X_umap_wnn"]          →  WNN UMAP
mdata.obsp["wnn_connectivities"]  →  WNN graph
mdata.obsp["wnn_distances"]       →  WNN distances
```

**Usage**

```python
from pipeline.modules.cite.cite_integration import run_mofa, run_totalvi

# MOFA+ (default, recommended for first run)
mdata, metrics = run_mofa(mdata, batch_key="batch", n_factors=15)
print(mdata.obsm["X_mofa"].shape)       # (21777, 15)
print(mdata.obsm["X_umap_mofa"].shape)  # (21777, 2)

# totalVI (better benchmark scores, slower)
mdata, metrics = run_totalvi(mdata, batch_key="batch", max_epochs=400)

# Quick test run
mdata, metrics = run_totalvi(mdata, batch_key="batch", max_epochs=10)
```

**Connects to**: `cite_integration_report.py` for the report

---

## 7. `cite_pipeline.py`

**Location**: `pipeline/modules/cite/cite_pipeline.py`

**What it does**

Orchestrates the full CITE-seq pipeline end-to-end in a single function call.
Calls all six Phase 4 modules in order, uses each step's output as input to the
next, and records per-step timing and metrics in a top-level provenance entry.

**Why it exists**

For use in notebooks or scripts where you want to run the full pipeline
programmatically without the CLI runner. Also useful for testing and
benchmarking. `run_cite_pipeline.py` (the CLI) calls this internally.

**Parameters (config dict)**

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `batch_key` | str | `"donor"` | obs column for Harmony + integration |
| `integration` | str | `"mofa"` | `"mofa"` or `"totalvi"` (case-insensitive) |
| `n_factors` | int | `15` | MOFA+ latent factors |
| `max_epochs` | int | `400` | totalVI training epochs |
| `annotation_map` | dict | None | Passed to `annotate_adt` |
| `filter_doublets` | bool | False | Passed to `detect_adt_doublets` |

**Input**

```
mdata  : MuData  →  must contain:
                    mdata["rna"]           = RNA with log1p in .X and raw in layers["counts"]
                    mdata["adt"]           = ADT with raw integer counts in .X
                    obs[batch_key]         = batch column on both modalities
```

**Output**

```
mdata  : MuData  →  fully processed MuData with all embeddings populated
                    mdata["adt"].layers["counts"]          = raw ADT counts
                    mdata["adt"].layers["adt_clr"]         = CLR-normalized ADT
                    mdata["adt"].obs["adt_doublet_score"]  = doublet scores
                    mdata["adt"].obsm["X_pca_adt"]         = ADT PCA
                    mdata["adt"].obsm["X_pca_harmony_adt"] = Harmony PCA
                    mdata["adt"].obsm["X_umap_adt"]        = ADT UMAP
                    mdata["adt"].obs["leiden"]             = Leiden clusters
                    mdata.obsm["X_mofa"]                   = MOFA+ factors (mofa)
                    mdata.obsm["X_umap_mofa"]              = MOFA+ UMAP (mofa)
                    mdata.uns["omicsage_cite_pipeline"]    = full provenance
```

**Usage**

```python
from pipeline.modules.cite.cite_pipeline import run_cite_pipeline

# Default config (MOFA+, batch_key="donor")
mdata_out = run_cite_pipeline(mdata)

# Custom config
mdata_out = run_cite_pipeline(mdata, config={
    "batch_key":    "batch",
    "integration":  "totalvi",
    "max_epochs":   50,
    "filter_doublets": True,
    "annotation_map": {"0": "CD4 T", "1": "B cell"},
})

# Access step timing
for step in mdata_out.uns["omicsage_cite_pipeline"]["steps"]:
    print(f"{step['step']:20s}  {step['elapsed_seconds']:.1f}s")
```

---

## 8. `run_cite_pipeline.py`

**Location**: repo root

**What it does**

Command-line runner for the CITE-seq pipeline. Mirrors the behaviour of
`run_pipeline.py` exactly — reads a YAML config, runs enabled steps in order,
checkpoints each step's output to disk, skips already-computed steps, and
generates an HTML report immediately after each step.

**Why it exists**

`cite_pipeline.py` is a Python function — it has no CLI, no checkpointing,
and no per-step reports. `run_cite_pipeline.py` adds all of that on top,
giving the analyst the same step-control workflow they have for the RNA pipeline.

**Step order**

```
normalize_adt → doublets → reduce_adt → harmony_adt → annotate_adt → integration
```

**Prerequisites**

The RNA pipeline must have been run first for the same dataset:

```bash
python run_pipeline.py --config config/runs/GSE194122.yaml
```

This produces the two files that the CITE pipeline needs:
- `data/processed/GSE194122/01_qc_adt.h5ad` — ADT input (steps 1–5)
- `data/processed/GSE194122/05_annotated.h5ad` — RNA input (step 6)

**CLI flags**

| Flag | Description |
|------|-------------|
| `--config PATH` | Required. Path to YAML config file |
| `--step NAME` | Run exactly one step |
| `--from-step NAME` | Start from this step (inclusive) |
| `--to-step NAME` | Stop at this step (inclusive) |
| `--force` | Re-run steps even if cached output exists |

`--from-step` and `--to-step` can be combined to run any contiguous range.
`--step` cannot be combined with `--from-step` / `--to-step`.

**Valid step names** (in order):
`normalize_adt` · `doublets` · `reduce_adt` · `harmony_adt` · `annotate_adt` · `integration`

**Checkpointing**

Each step writes its output to `data/processed/<dataset_id>/cite_NN_<step>.h5ad`.
Steps 1–5 write AnnData; step 6 (integration) writes MuData (`.h5mu`).
If the output file already exists the step is skipped and the cached path is
passed to the next step. Use `--force` to override.

**Reports**

An HTML report is generated immediately after each step completes, written to
`reports/<dataset_id>/cite_NN_<step>_report.html`. Report failures never crash
the pipeline — they are caught and printed as warnings.

**Validation**

Before any step runs, the runner checks that every required input file exists.
If a file is missing it exits with a clear message identifying which step needs
to be run first.

**Usage**

```bash
conda activate omicsage

# Full pipeline — all 6 steps
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml

# ADT steps only — skip integration
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --to-step annotate_adt

# Run from Harmony onward (normalize, doublets, reduce already done)
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --from-step harmony_adt

# Re-run annotation after filling in annotation_map in the config
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --step annotate_adt --force

# Switch integration from MOFA+ to totalVI — edit config, then:
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --step integration --force

# Run a range
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --from-step reduce_adt --to-step harmony_adt
```

**Interactive annotation workflow**

The standard workflow for a new dataset requires inspecting clusters before
applying labels. The runner supports this explicitly:

```bash
# 1. Run through annotation without an annotation_map
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --to-step annotate_adt

# 2. Open reports/GSE194122/cite_05_annotate_report.html
#    Inspect the cluster sizes table and UMAP
#    Identify each cluster from protein expression patterns

# 3. Add annotation_map to config under steps.annotate_adt.params:
#      annotation_map:
#        "0": "CD4 T"
#        "1": "CD8 T"
#        "2": "B cell"
#        "3": "NK"
#        "4": "Monocyte"

# 4. Re-run annotation with the map applied
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --step annotate_adt --force

# 5. Continue to integration
python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml \
    --step integration
```

---

## 9. `config/runs/GSE194122_cite.yaml`

**Location**: `config/runs/GSE194122_cite.yaml`

**What it does**

YAML configuration file for the CITE-seq pipeline on the GSE194122 BMMC
dataset. Specifies all step parameters, input paths, and output directories.
All fields have defaults — only change what differs from the defaults.

**Key paths**

| Key | Value | Description |
|-----|-------|-------------|
| `paths.adt_input` | `data/processed/GSE194122/01_qc_adt.h5ad` | ADT input — produced by RNA pipeline QC |
| `paths.rna_input` | `data/processed/GSE194122/05_annotated.h5ad` | Annotated RNA — produced by RNA pipeline annotate step |
| `paths.processed_dir` | `data/processed/GSE194122` | Where checkpoint h5ad files are written |
| `paths.reports_dir` | `reports/GSE194122` | Where HTML reports are written |

**Global batch key**

```yaml
cite:
  batch_key: batch
```

This is read by both `harmony_adt` and `integration` steps. Override per-step
under `steps.harmony_adt.params.batch_key` or `steps.integration.params.batch_key`
if needed.

**Common config patterns**

**No batch effects — skip Harmony:**
```yaml
harmony_adt:
  enabled: false
```
Note: `adt_reduce.py` still computes a neighbor graph. `adt_annotate.py` reads
`obsp["connectivities"]` from reduce, not harmony — so annotation still works.

**Skip integration entirely:**
```yaml
integration:
  enabled: false
```

**Use totalVI instead of MOFA+:**
```yaml
steps:
  integration:
    params:
      method: totalvi
      max_epochs: 400
```

**Apply annotation map:**
```yaml
steps:
  annotate_adt:
    params:
      annotation_map:
        "0": "CD4 T"
        "1": "CD8 T"
        "2": "B cell"
        "3": "NK"
        "4": "Monocyte"
        "5": "DC"
```

**Filter doublets (remove, not just flag):**
```yaml
steps:
  doublets:
    params:
      filter_doublets: true
```

**New dataset** — copy this file, update `dataset.id`, `paths.adt_input`,
`paths.rna_input`, and any parameters that differ. Zero Python required.

---

## 10. Report modules

All six report modules live in `reports/templates/cite/`. Each has a single
public function and follows the same pattern: takes the step's output AnnData
(or MuData for integration) plus the metrics dict, writes a self-contained
HTML file, and returns the absolute output path as a string.

All reports share the same visual style as the RNA pipeline reports:
dark-blue gradient header, stat cards with `#f0f2ff` background, white section
cards with subtle box-shadow, 130 dpi PNG figures embedded as base64.

### `cite_normalize_report.py`

**Function**: `run_cite_normalize_report(adt, metrics, report_path, dataset_name)`

| Figure | What it shows |
|--------|--------------|
| CLR violin | CLR value distribution per protein (up to 40 proteins shown) |
| Mean CLR bar | Mean CLR expression per protein, ranked highest to lowest |
| Raw vs CLR scatter | Per-cell raw total counts vs CLR sum — shows normalization effect |

---

### `cite_doublets_report.py`

**Function**: `run_cite_doublets_report(adt, metrics, report_path, dataset_name)`

| Figure | What it shows |
|--------|--------------|
| Score histogram | Distribution of `adt_doublet_score` with threshold line |
| Score by donor | Boxplot of doublet scores per donor/batch — spots problematic samples |

Also shows a summary table with pairs evaluated, pairs skipped, and whether
filtering was applied.

---

### `cite_reduce_report.py`

**Function**: `run_cite_reduce_report(adt, metrics, report_path, dataset_name)`

| Figure | What it shows |
|--------|--------------|
| ADT PCA | PC1 vs PC2 — plain scatter and coloured by batch side-by-side |
| ADT UMAP | Pre-Harmony UMAP coloured by batch — batch effects visible here |

The pre-Harmony UMAP is important as a baseline for comparing with the
post-Harmony UMAP in the next report.

---

### `cite_harmony_report.py`

**Function**: `run_cite_harmony_report(adt, metrics, report_path, dataset_name)`

| Figure | What it shows |
|--------|--------------|
| PCA before/after | Side-by-side: uncorrected `X_pca_adt` vs `X_pca_harmony_adt`, coloured by batch |
| Post-Harmony UMAP | `X_umap_adt` after correction — batch mixing should be visible |

The side-by-side PCA is the primary diagnostic: cells should be more intermixed
by batch in the right panel than the left.

---

### `cite_annotate_report.py`

**Function**: `run_cite_annotate_report(adt, metrics, report_path, dataset_name)`

| Figure | What it shows |
|--------|--------------|
| Cluster sizes | Bar chart of cells per Leiden cluster |
| UMAP by cluster | `X_umap_adt` coloured by Leiden ID — or by `adt_celltype` if annotation_map was applied |

Also shows a cluster table with cell counts and cell type labels (if annotated),
and a note reminding the analyst to fill in `annotation_map` when it has not
been provided.

---

### `cite_integration_report.py`

**Function**: `run_cite_integration_report(mdata, metrics, report_path, dataset_name)`

| Figure | What it shows |
|--------|--------------|
| UMAP by batch | Joint MOFA+/totalVI UMAP coloured by batch — should show mixing |
| UMAP by cell type | Same UMAP coloured by `cell_type_vote` from RNA (or `adt_celltype`) |
| UMAP by ADT cluster | Same UMAP coloured by ADT Leiden cluster |

Auto-detects which integration method was used from `mdata.obsm` keys
(`X_mofa` vs `X_totalVI`). All three panels use the same coordinate space
so correspondence between batch structure, cell type, and ADT cluster is
directly visible.

---

## 11. Checkpoint and report map

| Step | Checkpoint file | Report file |
|------|----------------|-------------|
| normalize_adt | `cite_01_normalized_adt.h5ad` | `cite_01_normalize_report.html` |
| doublets | `cite_02_doublets_adt.h5ad` | `cite_02_doublets_report.html` |
| reduce_adt | `cite_03_reduced_adt.h5ad` | `cite_03_reduce_report.html` |
| harmony_adt | `cite_04_harmony_adt.h5ad` | `cite_04_harmony_report.html` |
| annotate_adt | `cite_05_annotated_adt.h5ad` | `cite_05_annotate_report.html` |
| integration | `cite_06_integration.h5mu` | `cite_06_integration_report.html` |

All paths are under `paths.processed_dir` and `paths.reports_dir` from the config.

Steps 1–5 produce AnnData (`.h5ad`). Step 6 produces MuData (`.h5mu`) because
it joins both RNA and ADT into a single object.

---

## 12. Embedding key reference

All embedding keys follow a strict naming convention to prevent collisions when
RNA and ADT AnnData objects are joined into a MuData.

**On `mdata["rna"]` (from RNA pipeline)**

| Key | Contents | Written by |
|-----|----------|-----------|
| `obsm["X_pca"]` | RNA PCA | `reduce.py` |
| `obsm["X_pca_harmony"]` | RNA Harmony PCA | `harmony_correct.py` |
| `obsm["X_umap"]` | RNA UMAP | `reduce.py` |

**On `mdata["adt"]` (from CITE pipeline)**

| Key | Contents | Written by |
|-----|----------|-----------|
| `obsm["X_pca_adt"]` | ADT PCA (uncorrected) | `adt_reduce.py` |
| `obsm["X_pca_harmony_adt"]` | ADT Harmony PCA | `adt_harmony.py` |
| `obsm["X_umap_adt"]` | ADT UMAP (post-Harmony) | `adt_harmony.py` |

**On `mdata` (joint space)**

| Key | Contents | Written by |
|-----|----------|-----------|
| `obsm["X_mofa"]` | MOFA+ latent factors | `cite_integration.py` |
| `obsm["X_umap_mofa"]` | UMAP from MOFA+ | `cite_integration.py` |
| `obsm["X_totalVI"]` | totalVI latent representation | `cite_integration.py` |
| `obsm["X_umap_totalVI"]` | UMAP from totalVI | `cite_integration.py` |
| `obsm["X_umap_wnn"]` | WNN UMAP | deferred |

---

## 13. `cite_combined_report.py`

**Location**: `reports/templates/cite/cite_combined_report.py`

**What it does**

Reads all individual CITE-seq step HTML reports from `reports_dir` and
assembles them into a single self-contained **tabbed HTML file** with one
tab per pipeline step.

Only tabs for reports that actually exist on disk are shown. This means
stopping the pipeline at any step (e.g. `--to-step annotate_adt`) still
produces a valid combined report covering all completed steps — no errors,
no empty tabs.

**Why it exists**

After a full or partial pipeline run, the `reports/<dataset>/` directory
contains up to six separate HTML files. A biologist or collaborator has no
obvious entry point and must open each file individually. The combined report
provides a single file that is the complete CITE-seq analysis record, with a
progress bar showing how much of the pipeline has been run.

It is the CITE-seq equivalent of the RNA pipeline's `combined_report.py`
(which produces `00_combined_report.html`). The two combined reports coexist
in the same `reports/GSE194122/` directory without conflict — the CITE combined
report uses the `cite_00_` prefix to distinguish it.

**Design decisions**

- Zero changes to existing step report generators — the combiner reads their
  output after the fact rather than modifying how they write
- Extracts `<main>` content from each step report; falls back to `<body>`
  minus header/footer if no `<main>` tag is present
- Step-specific CSS is preserved by extracting and re-embedding `<style>` blocks
- No new dependencies — uses only Python stdlib (`re`, `datetime`, `pathlib`)
- Report failures are wrapped in `try/except` in `run_cite_pipeline.py` so
  a combiner failure never crashes the pipeline

**Tab registry** (in display order)

| Filename | Tab label | Icon |
|----------|-----------|------|
| `cite_01_normalize_report.html` | Normalize ADT | 🧪 |
| `cite_02_doublets_report.html` | Doublets | 🔍 |
| `cite_03_reduce_report.html` | Reduce ADT | 🔭 |
| `cite_04_harmony_report.html` | Harmony ADT | 🎵 |
| `cite_05_annotate_report.html` | Annotate ADT | 🏷️ |
| `cite_06_integration_report.html` | Integration | 🔗 |

**Parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reports_dir` | Path | required | Directory containing the CITE step HTML reports |
| `dataset_name` | str | `"CITE-seq Analysis"` | Shown in the report header and browser tab title |
| `output_path` | Path | `reports_dir / "cite_00_combined_report.html"` | Output path |

**Output**

```
cite_00_combined_report.html  →  single self-contained tabbed HTML file
                                  sorts to the top of the reports folder
                                  due to the "cite_00_" prefix
```

**Header features**

- Dataset name + generation timestamp
- Progress bar: `N of 6 pipeline steps complete` + percentage
- Step labels listed inline below the progress bar

**Navigation features**

- Click any tab to switch panels
- Left / right arrow keys navigate between tabs
- Tab button scrolls into view on narrow screens (horizontal scroll)

**Usage — called automatically by `run_cite_pipeline.py`**

The combined report is generated at the end of every `run_cite_pipeline.py`
run, after all active steps complete. It covers whichever steps ran — a
`--to-step annotate_adt` run produces a 5-tab combined report.

```python
from reports.templates.cite.cite_combined_report import generate_cite_combined_report

generate_cite_combined_report(
    reports_dir=Path("reports/GSE194122"),
    dataset_name="BMMC CITE-seq (NeurIPS 2021)",
)
# → reports/GSE194122/cite_00_combined_report.html
```

**Usage — standalone rebuild from existing reports**

If you want to regenerate the combined report without re-running the pipeline
(e.g. after manually editing a step report, or after copying reports to a new
location):

```bash
python -m reports.templates.cite.cite_combined_report \
    --reports-dir reports/GSE194122 \
    --dataset-name "BMMC CITE-seq (NeurIPS 2021)"

# Custom output path:
python -m reports.templates.cite.cite_combined_report \
    --reports-dir reports/GSE194122 \
    --dataset-name "BMMC CITE-seq (NeurIPS 2021)" \
    --output reports/GSE194122/cite_combined.html
```

**Coexistence with the RNA combined report**

Both combined reports live in the same `reports/GSE194122/` directory:

```
reports/GSE194122/
  00_combined_report.html        ←  RNA pipeline combined report
  cite_00_combined_report.html   ←  CITE-seq pipeline combined report  ← this file
  cite_01_normalize_report.html
  cite_02_doublets_report.html
  ...
```

The `cite_00_` prefix keeps the CITE combined report sorted immediately after
the RNA combined report when the directory is listed alphabetically.

**Connects to**: all six CITE-seq step report generators — reads their output
files, produces no `.h5ad` or `.h5mu` output of its own.
