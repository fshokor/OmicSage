# OmicSage — Annotation Module Extension Plan
> Status: DEFERRED — decisions locked, implementation pending
> Created: 2026-05-07
> Resume this plan in a dedicated session before DEG or clustering work expands further

---

## Decision: Drop rpy2 Permanently

rpy2 is not a viable dependency for OmicSage. Reasons:

- Requires system R 4.3+ installed separately from conda environment
- Requires `libtirpc-dev` system library (not standard on all Linux distros)
- Requires sudo access for R package installation
- Bioconductor packages (SingleR, celldex) download large references at runtime
- Breaks on HPC systems, Docker images, and collaborator machines without careful setup
- CI would require a full R installation — significant maintenance burden
- Violates the OmicSage principle: `pip install omicsage` should just work

**rpy2 will never appear in requirements.txt or requirements-ci.txt.**

The R methods (ScType, SingleR) will be reimplemented natively in Python.

---

## Method 3: ScType-py

### What ScType does
- Downloads ScTypeDB_full.xlsx (~150KB) from GitHub on first call
- Filters gene sets by tissue parameter (e.g. "Immune system", "Liver", "Brain")
- For each cluster, computes:
  `score = mean(positive marker expression) - mean(negative marker expression)`
- Assigns cluster → highest scoring cell type
- Writes `obs['cell_type_sctype']`

### Python implementation plan
```python
def _run_sctype(
    adata: AnnData,
    leiden_col: str,
    tissue: str = "Immune system",
    db_path: Path = Path("data/references/sctype/ScTypeDB_full.xlsx"),
) -> None:
    """
    Pure Python reimplementation of ScType scoring.
    Downloads ScTypeDB_full.xlsx on first call, caches locally.
    Writes obs['cell_type_sctype'].
    """
```

### Reference database
- File: ScTypeDB_full.xlsx
- URL: https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx
- Size: ~150 KB
- Format: Excel — columns: tissueType, cellName, geneSymbolmore (positive), shortGene (negative)
- Cache path: data/references/sctype/ScTypeDB_full.xlsx
- Download: once on first call via requests, skip if cached

### Dependencies
- `requests` — download DB (already in most envs)
- `openpyxl` — read xlsx (already installed for reports)
- `numpy` — scoring (already installed)
- No new dependencies needed

### Tissue values in ScTypeDB (confirmed)
"Immune system", "Liver", "Brain", "Lung", "Kidney", "Pancreas",
"Heart", "Intestine", "Skin", "Adrenal", "Bladder", "Breast",
"Cervix", "Eye", "Muscle", "Ovary", "Prostate", "Stomach", "Thyroid"

### Vote integration
- obs column: `obs['cell_type_sctype']`
- Already wired in `evidence_slots` in `_run_majority_vote()`:
  `("ScType", "cell_type_sctype", 1)`
- Double weight for parenchymal types already coded in `_PARENCHYMAL`
- No changes to vote logic needed — just implement the method

---

## Method 4: SingleR-py

### What SingleR does (algorithm)
1. For each pair of reference cell types → find pairwise marker genes
   (genes with high variance between those two types)
2. For each query cell → compute Spearman correlation vs each reference type
   using only the pairwise marker genes
3. Iterative refinement: keep top candidates, re-score with their specific markers
4. Pruning: low-confidence calls → "Unassigned"
5. Writes per-cell label → aggregate to cluster level

### Why not use the original R SingleR reference (HPCA)
- HumanPrimaryCellAtlasData() = ~200MB HDF5 download
- Too large to bundle, too slow to fetch at runtime
- Bulk RNA-seq from 2013 — lower resolution than modern scRNA-seq
- External Bioconductor dependency — fragile

### OmicSage compact reference strategy
Build `omicsage_reference.h5ad` — a compact pseudobulk reference:
- Source: well-annotated public scRNA-seq datasets (we already have two)
- Format: mean log-normalized expression per cell type (pseudobulk)
- Gene subset: top 200 highly variable genes per cell type
- Size: ~3-5 MB total (vs 200MB for HPCA)
- Expandable: any user with a labeled H5AD can contribute new cell types

### Reference datasets (build order)
| Dataset | Tissue | Cell types | Status |
|---------|--------|------------|--------|
| GSE194122 BMMC CITE-seq | Blood / Bone marrow | ~20 immune types | Already downloaded ✓ |
| GSE166635 HCC | Liver | ~11 types incl. hepatocyte | Already downloaded ✓ |
| TBD lung dataset | Lung | ~10 types | Future |
| TBD brain dataset | Brain | ~10 types | Future |
| TBD gut dataset | Gut epithelium | ~8 types | Future |

Coverage after first two datasets: blood + liver = ~50% of published scRNA-seq papers.

### Reference building script (to write)
```
scripts/build_singler_ref.py
    Input:  list of (h5ad_path, cell_type_col) tuples
    Output: data/references/singler/omicsage_ref_v1.h5ad
    Logic:  for each dataset → pseudobulk per cell type →
            select top 200 HVGs per type → concatenate → save
```

### Python implementation plan
```python
def _run_singler(
    adata: AnnData,
    leiden_col: str,
    ref_path: Optional[Path] = None,   # None → built-in lite ref
                                        # Path → user-supplied H5AD
                                        # "hpca" → download full HPCA (opt-in)
) -> None:
    """
    Lightweight Python reimplementation of SingleR correlation scoring.
    Uses compact pseudobulk reference by default.
    Writes obs['cell_type_singler'].
    """
```

### Dependencies
- `scipy.stats.spearmanr` — correlation scoring (already installed)
- `anndata` — load reference H5AD (already installed)
- `numpy` — matrix ops (already installed)
- No new dependencies needed

### Vote integration
- obs column: `obs['cell_type_singler']`
- Already wired in `evidence_slots` in `_run_majority_vote()`:
  `("SingleR_HPCA", "cell_type_singler", 1)`
- No changes to vote logic needed — just implement the method

---

## Updated annotate() signature (when implemented)

```python
adata_ann, ann_dict = annotate(
    adata_clustered,
    methods=["celltypist", "markers", "sctype", "singler", "vote"],
    leiden_col="leiden",
    celltypist_models=["Immune_All_High.pkl", "Immune_All_Low.pkl"],
    celltypist_models_dir=None,
    marker_sets=None,
    tissue="Immune system",        # ScType DB tissue filter
    singler_ref=None,              # None=built-in, Path=custom, "hpca"=download
    inplace=False,
)
```

New obs columns after full implementation:
```
obs['cell_type_sctype']    — ScType best label per cluster
obs['cell_type_singler']   — SingleR-py label per cluster
obs['cell_type_vote']      — upgraded to 4-way consensus
obs['cell_type_confidence']— fraction of 4 methods agreeing
```

---

## Tests to write (test_annotate.py additions)

All marked `@pytest.mark.skipif(os.getenv("OMICSAGE_CI") == "true", ...)`
since they require the reference files to be present.

```
test_sctype_column_exists()
test_sctype_labels_are_strings()
test_sctype_every_cluster_covered()
test_sctype_tissue_filter()         — "Liver" tissue gives hepatocyte types
test_singler_column_exists()
test_singler_labels_are_strings()
test_singler_every_cluster_covered()
test_singler_unassigned_for_low_confidence()
test_4way_vote_columns_exist()
test_4way_confidence_range()        — all values between 0.0 and 1.0
test_custom_singler_ref()           — user-supplied H5AD path works
```

---

## Session checklist (when resuming)

- [ ] Write `scripts/build_singler_ref.py`
- [ ] Run script on GSE194122 + GSE166635 → `omicsage_ref_v1.h5ad`
- [ ] Implement `_run_sctype()` in `pipeline/modules/qc/annotate.py`
- [ ] Implement `_run_singler()` in `pipeline/modules/qc/annotate.py`
- [ ] Add `tissue` and `singler_ref` params to `annotate()` signature
- [ ] Remove the "planned for Session B" warning block from `annotate()`
- [ ] Add 11 new tests to `tests/test_annotate.py`
- [ ] Update `reports/annotate_report.py` — add ScType + SingleR UMAP panels
- [ ] Update notebook Step 5 — methods list + sanity check cells
- [ ] Update `MODULE_DOCS.md`

---

## Notes
- ScType can be done in one session (~3 hours) — simpler algorithm, DB already known
- SingleR-py needs two sessions: one for `build_singler_ref.py`, one for the scorer
- Do ScType first, SingleR second
- Resume AFTER DEG module is complete — annotation quality improves
  once we can validate labels against DEG marker genes
