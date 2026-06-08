# NEXT_SESSION.md — Phase 7 Session 2
# Spatial Transcriptomics: Normalization, HVG, PCA, Spatial Neighbours
> Updated: June 2026
> Previous session: Phase 7 Session 1 — spatial_ingest.py + spatial_qc.py COMPLETE

---

## Session Context

**Phase:** 7 — Spatial Transcriptomics
**Session number:** 2 of 4
**Last thing completed:** Phase 7 Session 1 —
  - `spatial_ingest.py`: unified dispatcher for visium / h5ad / benchmark;
    stubs for xenium / merfish / codex / visium_hd with NotImplementedError;
    auto-detection fingerprints; `list_supported_types()`
  - `spatial_qc.py`: spot filtering on total_counts, n_genes, MT%;
    `qc_pass` mask; provenance in `uns["omicsage_spatial_qc"]`
  - `spatial_qc_report.py`: self-contained HTML with violin plots,
    spatial scatter (total_counts + MT%), threshold bar chart
  - `run_spatial_pipeline.py`: stub runner (Sessions 1 steps only)
  - `config/runs/visium_lymphnode.yaml`
  - `tests/test_spatial_ingest_qc.py`

**Key architectural decision made this session:**
  - `spatial_type` (not `source_type`) is the config/API key for technology
  - `_LOADER_REGISTRY` + `_AUTO_FINGERPRINTS` pattern for extensibility
  - All downstream modules receive a standard AnnData — technology is
    transparent to everything after ingest


---

## Spatial Pipeline Architecture (reference for all Phase 7 sessions)

```
┌─────────────────────────────────────────────────────────────────┐
│  INGEST LAYER — spatial_ingest.py (one file, technology-aware)  │
│                                                                 │
│  spatial_type="visium"    → sq.read.visium()        ✅ done     │
│  spatial_type="h5ad"      → sc.read_h5ad()          ✅ done     │
│  spatial_type="benchmark" → sq.datasets.*           ✅ done     │
│  spatial_type="visium_hd" → sq.read.visium_hd()     🔲 planned  │
│  spatial_type="xenium"    → sq.read.xenium()        🔲 planned  │
│  spatial_type="merfish"   → sq.read.vizgen()        🔲 planned  │
│  spatial_type="codex"     → custom CSV loader       🔲 planned  │
│                                                                 │
│  Auto-detection from directory fingerprints.                    │
│  Config key: spatial.spatial_type (default: "auto")            │
└────────────────────────────┬────────────────────────────────────┘
                             │ standard AnnData
                             │ obsm["spatial"], uns["spatial"]
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  QC LAYER                                                       │
│                                                                 │
│  spatial_qc.py     → counts, n_genes, MT% filtering  ✅ done    │
│  (future) spatial_qc_imaging.py → cell area, etc.    🔲 future  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  REDUCTION LAYER — shared across all technologies               │
│                                                                 │
│  spatial_reduce.py → normalize, HVG, PCA,                      │
│                       sq.gr.spatial_neighbors    ← SESSION 2    │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  CLUSTERING LAYER — shared across all technologies              │
│                                                                 │
│  spatial_cluster.py → Leiden + spatially variable genes         │
│                        (Moran's I)               ← SESSION 3    │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  DECONVOLUTION LAYER — Visium / Visium HD only                  │
│  (skipped for single-cell resolution technologies)              │
│                                                                 │
│  spatial_deconvolve.py → cell2location + report  ← SESSION 4    │
└─────────────────────────────────────────────────────────────────┘
```

**Key design principle:** `spatial_ingest.py` normalises all technologies
into a standard AnnData with `obsm["spatial"]`. Everything from
`spatial_reduce.py` downward is technology-agnostic — it never checks
`spatial_type`. The technology label flows only through provenance
(`uns["omicsage_spatial_ingest"]["spatial_type"]`), used for report
titles and to conditionally skip deconvolution for single-cell-resolution
data.

**Current test count:** 1105+ prior + new spatial tests = ~1160+
(run baseline to confirm exact count before starting)

---

## Baseline Verification

```bash
conda activate omicsage
cd ~/OmicSage
python -m pytest tests/ -q --tb=short 2>&1 | tail -5
```

Expected: 1160+ passed, 0 errors.

Also verify the spatial modules import cleanly:
```bash
python -c "
from pipeline.modules.spatial.spatial_ingest import spatial_ingest, list_supported_types
from pipeline.modules.spatial.spatial_qc import spatial_qc
print('imports OK')
print('supported types:', list_supported_types())
"
```

---

## Today's Goal — `spatial_reduce.py`

Normalization, HVG selection, PCA, and spatial neighbours graph for
Visium spots. This is the bridge between raw QC-filtered counts and
the clustering step in Session 3.

Input:  AnnData from `spatial_qc` (has `obsm["spatial"]`, `qc_pass` column)
Output: AnnData with `X` normalized, `obsm["X_pca"]`, spatial graph in
        `obsp["spatial_connectivities"]`, provenance in
        `uns["omicsage_spatial_reduce"]`

---

## Planned steps this phase (updated)

| Session | Deliverable | Status |
|---|---|---|
| 1 | `spatial_ingest.py` + `spatial_qc.py` + report + tests | ✅ DONE |
| **2 (this)** | `spatial_reduce.py` — normalize, HVG, PCA, spatial neighbours | ← |
| 3 | `spatial_cluster.py` — Leiden clustering + spatially variable genes (Moran's I) | |
| 4 | `spatial_deconvolve.py` — cell2location + combined HTML report | |

---

## Function signature (proposed)

```python
def spatial_reduce(
    adata: AnnData,
    n_top_genes: int = 3000,          # HVGs to select
    n_comps: int = 50,                # PCA components
    n_neighbors: int = 6,             # spatial graph neighbours (Visium hex grid = 6)
    coord_type: str = "grid",         # "grid" (Visium) or "generic"
    normalize_total: bool = True,     # sc.pp.normalize_total
    target_sum: float = 1e4,          # normalization target
    log1p: bool = True,               # sc.pp.log1p after normalization
    flavor: str = "seurat",           # HVG flavor: "seurat" | "cell_ranger"
    inplace: bool = False,
) -> tuple[AnnData, dict]:
```

---

## Key implementation decisions to confirm before coding

1. **`coord_type="grid"` vs `"generic"`** — Visium spots sit on a hexagonal
   grid so `sq.gr.spatial_neighbors(coord_type="grid")` gives the correct
   6-neighbour connectivity. Confirm this is the right call for standard
   Visium (not Visium HD or other technologies).

2. **HVG flavor** — `flavor="seurat"` calls `expm1` internally which can
   overflow on TF-IDF ATAC data. Visium RNA data is not TF-IDF so seurat
   is safe here. Confirm.

3. **Raw counts** — normalization must be applied to raw counts, not
   already-normalized data. After `spatial_qc`, `adata.X` is still raw
   (we never normalize in QC). Verify this assumption holds for benchmark
   data loaded via `sq.datasets.visium_hne_adata()` — the benchmark
   dataset is pre-processed (already log-normalized). Handle both cases:
   check `adata.uns.get("omicsage_spatial_ingest", {}).get("spatial_type")`
   and skip normalization if data is already normalized.
   SAFER APPROACH: always store raw in `layers["counts"]` at ingest and
   normalize from there, resetting `adata.X`.

4. **`n_neighbors=6`** — standard for Visium hex grid. Configurable for
   non-Visium data.

5. **PCA on HVGs only** — subset `adata[:, adata.var.highly_variable]`
   before PCA, same as scRNA pipeline.

---

## Pre-implementation research required

Before writing any code, fetch and read:

1. **squidpy spatial_neighbors API** — confirm `coord_type`, `n_neighs`
   parameters and what keys it writes to `obsp` and `uns`:
   - https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html

2. **sc-best-practices spatial normalization** — confirm normalize_total +
   log1p is recommended (not SCTransform) for Python/squidpy workflow:
   - https://www.sc-best-practices.org/spatial/spatially_variable_genes.html
   - https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html

3. **HVG selection on Visium** — confirm `sc.pp.highly_variable_genes`
   with `flavor="seurat"` is appropriate and n_top_genes=3000 is
   reasonable for Visium (~33k genes):
   - https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html

---

## Files to create this session

```
pipeline/modules/spatial/spatial_reduce.py      ← new
reports/templates/spatial/spatial_reduce_report.py  ← new (embed PCA variance + HVG plots)
tests/test_spatial_reduce.py                    ← new
```

Update:
```
run_spatial_pipeline.py                         ← add spatial_reduce step
```

---

## Expected AnnData state after spatial_reduce

```
adata.X                          normalized + log1p counts
adata.layers["counts"]           raw integer counts (preserved from ingest)
adata.var["highly_variable"]     bool mask of selected HVGs
adata.var["means"]               gene means (from HVG selection)
adata.var["dispersions_norm"]    normalized dispersions
adata.obsm["X_pca"]              PCA embedding (n_obs × n_comps)
adata.uns["pca"]["variance_ratio"] explained variance per PC
adata.obsp["spatial_connectivities"] spatial graph adjacency (from sq.gr.spatial_neighbors)
adata.obsp["spatial_distances"]  spatial distances
adata.uns["spatial_neighbors"]   spatial graph metadata
adata.uns["omicsage_spatial_reduce"] provenance dict
```

---

## Report content for spatial_reduce_report.py

- HVG scatter plot: mean vs dispersion, HVGs highlighted
  (`sc.pl.highly_variable_genes`)
- PCA variance explained: elbow plot (cumulative variance ratio per PC)
- PCA scatter: obs coloured by `total_counts` (spot quality check)
- Spatial connectivity: number of neighbours per spot (should be ~6 for Visium)
- Summary table: n_hvgs selected, n_pcs, spatial graph stats

---

## Known issues / watch-outs from Session 1

- Benchmark dataset (`sq.datasets.visium_hne_adata()`) is **pre-processed**
  (already normalized + log1p + clustered). Do NOT re-normalize it.
  Detect by checking `adata.uns.get("omicsage_spatial_ingest", {})
  .get("spatial_type") == "benchmark"` and skip normalization,
  OR better: always work from `layers["counts"]` if present.

- `sc.read_visium` is deprecated since scanpy 1.11 — always use
  `sq.read.visium()`. Already handled in `spatial_ingest.py`.

- `flavor="seurat"` for HVG calls `expm1` — safe for RNA, would overflow
  on TF-IDF ATAC. Not an issue here (spatial RNA only).

---

## Phase 7 full reference list

### Primary references (read before each session)

| Session | URL | What to read |
|---|---|---|
| 1 ✅ | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html | Full Visium tutorial — loading, QC, analysis |
| 1 ✅ | https://squidpy.readthedocs.io/en/stable/api/squidpy.read.visium.html | `sq.read.visium()` API |
| 1 ✅ | https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_visium.html | `sc.read_visium()` deprecation notice |
| **2** | https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html | `sq.gr.spatial_neighbors()` API — coord_type, n_neighs |
| **2** | https://www.sc-best-practices.org/spatial/spatially_variable_genes.html | Normalization + HVG + PCA for spatial data |
| **2** | https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html | HVG selection parameters |
| 3 | https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_autocorr.html | Moran's I for spatially variable genes |
| 3 | https://www.sc-best-practices.org/spatial/spatially_variable_genes.html | SVG methods comparison (Moran's I vs SpatialDE) |
| 4 | https://www.sc-best-practices.org/spatial/deconvolution.html | cell2location tutorial (full) |
| 4 | https://cell2location.readthedocs.io/en/latest/ | cell2location API docs |
| 4 | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html | Neighbourhood enrichment + co-occurrence (downstream) |

### Secondary references (background reading)

- **Squidpy paper**: Palla et al., Nature Methods 2022
  https://doi.org/10.1038/s41592-021-01358-2
  — Framework paper; cites all the spatial statistics methods used in Sessions 2-3

- **cell2location paper**: Kleshchevnikov et al., Nature Biotechnology 2022
  https://doi.org/10.1038/s41587-021-01139-4
  — Read before Session 4; explains N_cells_per_location hyperparameter

- **Deconvolution benchmark**: Li et al., Nature Methods 2022
  https://doi.org/10.1038/s41592-022-01480-9
  — Shows cell2location performs best on simulated data; RCTD best on real data
  — Useful for DECISIONS.md entry on why we chose cell2location

- **10x Visium Space Ranger output format**:
  https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview
  — Reference for understanding what files `sq.read.visium()` reads

- **OSTA book** (R-focused but conceptually solid):
  https://lmweber.org/OSTA-book/
  — Good background on spatially variable genes and clustering concepts

---

## DECISIONS.md entries to add this session

After Session 2, add to `DECISIONS.md`:

```markdown
## Phase 7 — Spatial Transcriptomics

### Normalization: normalize_total + log1p (not SCTransform)
- SCTransform is R-only (Seurat). Python/scanpy ecosystem uses
  normalize_total + log1p as the standard.
- sc-best-practices spatial chapter confirms this approach.

### HVG: n_top_genes=3000, flavor="seurat"
- Visium has ~33k genes; 3000 HVGs is standard for spatial.
- flavor="seurat" is safe for RNA data (no TF-IDF overflow risk unlike ATAC).

### Spatial neighbours: coord_type="grid", n_neighs=6
- Visium spots sit on a hexagonal grid — 6 neighbours is the correct
  connectivity for standard Visium.
- For non-Visium data (Xenium, MERFISH etc.) use coord_type="generic".

### Deconvolution: cell2location (not RCTD)
- sc-best-practices recommends cell2location; best on simulated benchmarks.
- Requires GPU for large datasets (>5k spots); CPU fallback available.
- Requires paired scRNA-seq reference — documented as user requirement.
- RCTD is simpler and CPU-only; noted as alternative in config comments.
```
