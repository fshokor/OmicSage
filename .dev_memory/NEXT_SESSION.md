# NEXT_SESSION.md — Phase 7 Session 1
# Spatial Transcriptomics: Visium data ingestion + QC
> Updated: June 2026
> Previous session: Phase 6 Session 1 — GRN inference COMPLETE (with known limitations)

---

## Session Context

**Phase:** 7 — Spatial Transcriptomics
**Session number:** 1 of ~4
**Last thing completed:** Phase 6 — multiome_grn.py running end-to-end,
  1176 ATAC TFs, 9695 edges, report renders. Known limitation: ATAC scores
  are flat (no coordinate-level motif matching). Documented in DECISIONS.md.

**Current test count:** 51 passed (test_multiome_grn.py) + 1054+ prior = ~1105+

**GRN known limitations (see DECISIONS.md):**
- ATAC scores all 1.0 — no cisTarget DB, no coordinate-level motif scan
- RNA TFs = 0 — pyscenic AUCell scoring not returning results
- Both deferred; acceptable for portfolio demo

---

## Baseline Verification

```bash
conda activate omicsage
cd ~/OmicSage
python -m pytest tests/ -q --tb=short 2>&1 | tail -5
```

Expected: 1105+ passed, 0 errors.

---

## Today's Goal — `spatial_qc.py`

Visium spatial transcriptomics QC module.

Input: 10x Visium Space Ranger output directory (or pre-built AnnData)
Output: AnnData with spatial coordinates, QC metrics, filtered spots

---

## Benchmark dataset

Use a publicly available 10x Visium human PBMC or FFPE tissue sample.
Recommended: 10x Genomics public dataset — Human Lymph Node
- Download: https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0
- Files needed: filtered_feature_bc_matrix/, spatial/

Alternatively use squidpy's built-in Visium dataset for faster setup:
```python
import squidpy as sq
adata = sq.datasets.visium_hne_adata()
```

---

## Planned steps this phase

| Session | Deliverable |
|---|---|
| **1 (this)** | `spatial_qc.py` — spot QC, MT%, spatial coordinate loading |
| 2 | `spatial_reduce.py` — HVG, PCA, spatial neighbors (squidpy) |
| 3 | `spatial_cluster.py` — Leiden + spatially variable genes |
| 4 | `spatial_deconvolve.py` — cell type deconvolution (RCTD/cell2location) + report |

---

## Function signature (proposed)

```python
def spatial_qc(
    adata: AnnData,                     # raw Visium AnnData (squidpy.read.visium or sc.read_visium)
    min_counts: int = 500,
    max_counts: int = 100000,
    min_genes: int = 200,
    max_genes: int = 10000,
    max_mt_pct: float = 20.0,
    mt_prefix: str = "MT-",
    filter_spots: bool = True,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
```

---

## Pre-implementation research required

Before writing any code:

1. **squidpy API** — confirm `sq.read.visium()` vs `sc.read_visium()` for loading
   - Reference: https://squidpy.readthedocs.io/en/stable/
   - sc-best-practices spatial chapter: https://www.sc-best-practices.org/spatial/spatially_variable_genes.html

2. **Spatial coordinates** — confirm obs keys: `obsm["spatial"]` or `obs["x_array"]`?

3. **QC metrics** — confirm `sc.pp.calculate_qc_metrics` works on Visium AnnData
   with spatial coords intact

4. **squidpy availability** — check if already in omicsage env:
   ```bash
   python -c "import squidpy; print(squidpy.__version__)"
   ```

---

## Key architectural decisions to make

- Use `squidpy.read.visium()` or `scanpy.read_visium()`? (squidpy wraps scanpy's)
- Store spatial coords in `obsm["spatial"]` (squidpy default) — confirm this
- QC plot: spatial scatter coloured by total_counts and pct_mt (squidpy plot)
- Report: HTML with spatial QC plots embedded as base64

---

## Files to create this session

```
pipeline/modules/spatial/spatial_qc.py
reports/templates/spatial/spatial_qc_report.py
tests/test_spatial_qc.py
run_spatial_pipeline.py               (stub — add steps as phase progresses)
config/runs/visium_lymphnode.yaml      (or squidpy built-in dataset)
```

---

## Reference URLs

- sc-best-practices spatial: https://www.sc-best-practices.org/spatial/spatially_variable_genes.html
- squidpy docs: https://squidpy.readthedocs.io/en/stable/
- squidpy spatial QC tutorial: https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html
- 10x Visium human lymph node dataset: https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0
- OSTA book (spatial transcriptomics): https://lmweber.org/OSTA-book/
