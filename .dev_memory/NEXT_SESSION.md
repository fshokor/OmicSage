# OmicSage — Next Session
> Written: 2026-06-11
> Phase: 7 extension — Spatial Imputation + New Format Support

---

## Session Context

**Last thing completed this session:**
- Fixed spatial combined report image presentation: lightbox click-to-expand,
  `max-width: 100%` on `.fig-wrap`, `max-width: 1400px` tab panel, `cursor: zoom-in`
  on all figures. Fixed `_stitch_panels` panoramic image bug (5582px wide → capped
  at 1800px via `TARGET_W` scale factor). All 6 spatial report files updated.
- Written `ai/manual_review/SPATIAL_MASTER_PROMPT.md` — 9-task spatial review
  prompt following the same structure as RNA/CITE/Multiome master prompts.

**Files modified this session:**
```
reports/templates/spatial/spatial_combined_report.py   ← lightbox + CSS fixes
reports/templates/spatial/spatial_qc_report.py         ← CSS max-width fix
reports/templates/spatial/spatial_reduce_report.py     ← CSS max-width fix
reports/templates/spatial/spatial_cluster_report.py    ← CSS + stitch figsize fix
reports/templates/spatial/spatial_deconvolve_report.py ← CSS + stitch figsize fix
reports/templates/spatial/spatial_downstream_report.py ← CSS max-width fix
ai/manual_review/SPATIAL_MASTER_PROMPT.md              ← NEW
```

**Verify last session works:**
```bash
conda activate omicsage
python -m pytest tests/test_spatial_*.py -q --tb=short
```
Confirm baseline test count before writing any new code.

---

## Two-Session Plan

This is a two-session extension of Phase 7. Do ONE session at a time.

---

## Session A (THIS SESSION) — Spatial Imputation

### Goal
Build `spatial_impute.py` — impute full transcriptome onto Visium spots using
a paired scRNA-seq reference. This is the missing chapter from sc-best-practices.

### Deliverables (all three together, in one session)
1. `pipeline/modules/spatial/spatial_impute.py`
2. `reports/templates/spatial/spatial_impute_report.py`
3. `tests/test_spatial_impute.py`
4. Update `run_spatial_pipeline.py` — add `"impute"` step
5. Update `reports/spatial_combined_report.py` — add Impute tab to TAB_REGISTRY
6. Update `config/runs/kuppe_heart.yaml` — add impute block

### Method Stack
- **Primary: Tangram** (`tangram-sc` on PyPI) — optimal transport sc→spatial mapping
- **Opt-in: gimVI** via `scvi-tools` — deep generative model, slower, higher memory
- Install: `pip install tangram-sc` (add to `environment.yml` + `requirements-ci.txt`)

### Module API
```python
def spatial_impute(
    adata_spatial: AnnData,      # Visium AnnData — post cluster checkpoint
    adata_sc: AnnData,           # paired scRNA-seq reference (annotated)
    method: str = "tangram",     # "tangram" | "gimvi"
    cell_type_key: str = "cell_type",
    n_top_genes: int = 2000,
    device: str = "cpu",
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    # Output:
    #   adata.obsm["imputed_expression"]  — DataFrame (spots × n_top_genes)
    #   adata.uns["omicsage_impute"]      — provenance dict
```

### Config block to add
```yaml
spatial:
  impute:
    enabled: true
    method: tangram            # tangram | gimvi
    sc_reference_path: ""      # path to paired scRNA-seq .h5ad
    n_top_genes: 2000
    cell_type_key: "cell_type"
    device: cpu
```

### Report sections (`spatial_impute_report.py`)
1. **Run Summary** — stat cards: n_genes_imputed, method, mapping_score (mean),
   n_spots, sc_reference (filename)
2. **Mapping Score Distribution** — histogram of per-spot Tangram mapping scores
   (quality metric: spots with score < 0.1 are poorly mapped)
3. **Top Imputed Genes on H&E** — spatial scatter of top 5 imputed genes by
   variance (the visual payoff — genes not in the original HVG set)
4. **Imputation Validation** — scatter plot: measured vs imputed expression for
   the genes present in both datasets (Spearman r in title)

### Runner step to add in `run_spatial_pipeline.py`
- Step name: `"impute"`
- Predecessor: `"cluster"` (does not require deconvolution)
- Skip condition: `config.get("spatial", {}).get("impute", {}).get("enabled", False) is False`
  OR `sc_reference_path` is empty/null
- Input checkpoint: cluster h5ad
- Output checkpoint: impute h5ad

### Test strategy (`test_spatial_impute.py`)
- Mock tangram import for unit tests (avoid requiring GPU/heavy install in CI)
- Test with tiny synthetic spatial AnnData (20 spots × 50 genes) +
  tiny sc reference (100 cells × 50 genes, 3 cell types)
- Test: output contract — `obsm["imputed_expression"]` present, correct shape
- Test: provenance in `uns["omicsage_impute"]`
- Test: `inplace=False` does not mutate input
- Test: `method="gimvi"` raises `ImportError` gracefully when scvi-tools absent
- Test: skips cleanly when `sc_reference_path` is null in config
- Test: report renders without error on mock data
- **Do NOT require real Tangram install for CI** — use `pytest.importorskip`
  pattern or mock

### Benchmark dataset note
The Kuppe heart dataset has a paired snRNA-seq reference:
`data/benchmark/kuppe_heart/snrna_ref.h5ad` (or equivalent).
Use this for manual end-to-end validation after tests pass.
If the paired reference is not yet downloaded, download it before starting:
- Kuppe et al. 2022 supplementary data — snRNA-seq processed h5ad
- GEO accession: GSE183852

---

## Session B (NEXT SESSION AFTER A) — New Format Support

### Goal
Activate the Xenium and Visium HD stubs already in `spatial_ingest.py`.
MERFISH and CODEX remain as `NotImplementedError` stubs.

### Why these two formats
- Visium HD: same 10x ecosystem, `sq.read.visium_hd()` reader exists, natural
  upgrade from Visium, public benchmark datasets available
- Xenium: most common new imaging-based platform, `sq.read.xenium()` reader
  exists in squidpy, cell-level (not spot-level) data model

### What changes per format

| Format | Key difference | Pipeline impact |
|--------|---------------|-----------------|
| Visium HD | No H&E by default; binned 2/8/16µm | `bin_size` config param; QC min_counts threshold lower |
| Xenium | Cell-level not spot-level; segmentation mask; transcripts.parquet | Cell-level QC; no `in_tissue` flag; different n_genes range |

**Everything downstream of ingest is format-agnostic** — QC, reduce, cluster,
deconvolve, impute all operate on the AnnData contract and need no changes.

### Files to change (Session B only)
```
pipeline/modules/spatial/spatial_ingest.py   ← implement _load_visium_hd + _load_xenium
tests/test_spatial_ingest_qc.py              ← stub tests → real tests
config/runs/                                 ← new yaml per format
SPATIAL_MODULE_DOCS.md                       ← format section updated
```

### Config additions for Session B
```yaml
spatial:
  source: "/path/to/data"
  spatial_type: "auto"     # auto | visium | visium_hd | xenium | h5ad | benchmark
  bin_size: 16             # Visium HD only: 2 | 8 | 16 µm
```

### Benchmark datasets for Session B (download before starting)
- **Visium HD:** 10x Genomics public dataset — Human Colorectal Cancer
  https://www.10xgenomics.com/datasets (search "Visium HD")
- **Xenium:** 10x Genomics public dataset — Human Breast Cancer
  https://www.10xgenomics.com/datasets (search "Xenium")
Both are free to download without account.

### Pre-Session B research step (mandatory)
Before writing any code for Xenium or Visium HD loaders, fetch the squidpy
source for both readers:
- `sq.read.visium_hd` — https://github.com/scverse/squidpy/blob/main/src/squidpy/read/_read.py
- `sq.read.xenium` — same file
Read the actual parameter names. Do NOT guess. This is the same rule that
caused bugs in Phase 7 Session 5 with `img_res_key` vs `img_key`.

---

## Known Issues / Watch Out For

- Tangram requires `torch` — confirm it is in `environment.yml` (it should be,
  scvi-tools depends on it). If not, add before starting.
- Tangram's `tg.pp_adatas()` mutates both input AnnData objects in place —
  work on copies inside the module.
- The Kuppe heart snRNA-seq reference may use `cell_type` or `cell_type_vote`
  as the annotation column — check before passing `cell_type_key`.
- `obsm["imputed_expression"]` is a DataFrame (genes as columns) — serialize
  to checkpoint as JSON string (same pattern as other DataFrame obsm keys in
  the spatial pipeline) to avoid h5py mixed-type errors.

---

## Memory Files to Update at End of Session A
```
.dev_memory/NEXT_SESSION.md    ← replace with Session B content
.dev_memory/CURRENT_STATUS.md  ← add impute to spatial pipeline complete list
.dev_memory/PROGRESS.md        ← tick spatial_impute.py
SPATIAL_MODULE_DOCS.md         ← add imputation section
```
