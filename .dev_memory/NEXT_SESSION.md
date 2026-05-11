## Session Context
Date: 2026-05-11
Phase: 1 — Core scRNA Pipeline
Last thing completed: Clustering on harmony-corrected embedding —
                      cluster.py extended (neighbors_key, cluster_key, compute_ari()),
                      test_cluster.py updated (16 tests), harmony_report.py UMAP fix,
                      notebook Steps 8 + 9 added, all four memory files updated.
File last worked on: .dev_memory/ (memory file updates)

## Today's Goal
Pseudobulk DEG — aggregate counts per (cell_type, donor) and run
pydeseq2 one-vs-rest per cell type.
New file: pipeline/modules/downstream/pseudobulk_deg.py
Output stored in the same format as deg.py so deg_report.py works unchanged.
ONE goal — do not start scVI, ADT/ATAC modules, or Quarto reports this session.

## Step 1 — Verify all prior tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_deg.py tests/test_gsea.py tests/test_harmony.py tests/test_cluster.py -v
# Expected: 11 + 8 + 13 + 16 = 48 passed
```

## Step 2 — Install pydeseq2 if not already present
```bash
conda activate omicsage
pip install pydeseq2
python -c "import pydeseq2; print(pydeseq2.__version__)"
```

## Step 3 — Input file for real-data run
```bash
# Input:  data/processed/GSE194122_cite_annotated.h5ad
#         (has obs['cell_type_vote'], obs['batch'], layers['counts'])
# Output: data/processed/GSE194122_cite_pseudobulk_deg.h5ad
```

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- OldFormatWarning from GSE194122 — harmless
- Notebook must be opened from OmicSage root
- seurat_v3 HVG flavor needs >=500 cells
- numpy bool: use bool() cast before `is True` checks in tests
- float dict keys must be stringified before storing in uns
- obs['cell_type_vote'] is the consensus column (not obs['cell_type'])
- obs['cell_type_groundtruth'] holds the preserved publication ground-truth
- rpy2 NOT used — do not attempt R integration
- gseapy.enrichr returns NO Overlap column in v1.1.3 — derived from Genes column
- deg.py requires rankby_abs=True in rank_genes_groups
- deg.py n_genes default is 500 (not 200)
- gsea.py direction="both" stores results under "{group}__up" / "{group}__down" keys
- harmonypy ho.Z_corr is (n_cells, n_pcs) — do NOT transpose
- harmony_correct() stores X_umap_harmony (not X_umap) and preserves
  original UMAP as X_umap_precorrection
- obs[batch_key] must be cast .astype(str) before pandas .map() or
  .value_counts() — Categorical + MultiIndex causes NotImplementedError
- reports/figures/ directory must be created before plt.savefig()
  (use Path("reports/figures").mkdir(parents=True, exist_ok=True))
- sc.pl.umap() reads obsm['X_umap'] by default — temporarily set
  adata.obsm["X_umap"] = adata.obsm["X_umap_harmony"] before calling it
- cluster.py _res_key() now takes (cluster_key, res) — not just (res)
- harmony_report.py umap_key default is now 'X_umap_harmony' (not 'X_umap')
- pseudobulk requires layers['counts'] (raw integers) — NOT layers['logcounts']
- pydeseq2 needs min_samples >= 2 per group — filter cell types with only 1 donor

## Deferred Work (do not start this session)
- scVI batch correction -> deferred to Phase 6 (MultiVI for multiome integration)
- ScType-py + SingleR-py annotation -> docs/ANNOTATION_PLAN.md
- ADT and ATAC QC modules
- Quarto report templates (Phase 2)

## Files Modified Last Session
### Updated files
- pipeline/modules/clustering/cluster.py   <- neighbors_key, cluster_key, compute_ari()
- tests/test_cluster.py                    <- 16 tests (was 12), _res_key signature fixed
- tests/test_harmony.py                    <- 13 tests (was 12), adata_with_umap fixture
- reports/harmony_report.py                <- before/after now both UMAPs
- notebooks/phase1_qc.ipynb                <- Steps 8 + 9 added (19 new cells)

### Memory files
- .dev_memory/MODULE_DOCS.md  <- Modules 19-20 added, data flow updated
- .dev_memory/CURRENT_STATUS.md <- Updated
- .dev_memory/NEXT_SESSION.md <- Updated (this file)
- .dev_memory/PROGRESS.md <- Updated

## Verify Last Session Works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_cluster.py -v
# Expected: 16 passed
```

## Relevant Context — GSE194122 Data Chain
- data/processed/GSE194122_cite_annotated.h5ad        <- input for pseudobulk_deg()
- data/processed/GSE194122_cite_harmony.h5ad           <- harmony output
- data/processed/GSE194122_cite_harmony_clustered.h5ad <- harmony + clustering output

## Pseudobulk DEG — design notes
- Aggregate layers['counts'] per (cell_type_vote, batch) -> bulk-like matrix
- Filter: min_cells=10 per group-donor combo, min_samples=3 donors per cell type
- Run pydeseq2 one-vs-rest per cell type (same grouping as deg.py)
- Return format identical to deg.py deg_dict so deg_report.py works unchanged
- Store provenance in uns['omicsage_pseudobulk_deg']
- Edge cases: cell types with < min_samples donors -> skip with warning, not crash

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, rpy2 3.6.7 (installed but NOT used)
New package needed: pydeseq2 (pip install pydeseq2)
