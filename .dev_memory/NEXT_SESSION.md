## Session Context
Date: 2026-05-11
Phase: 1 — Core scRNA Pipeline
Last thing completed: Pseudobulk DEG —
                      pseudobulk_deg.py, test_pseudobulk_deg.py (14 tests),
                      pseudobulk_deg_report.py, notebook Step 10 added,
                      all four memory files updated.
File last worked on: .dev_memory/ (memory file updates)

## Today's Goal
[DECIDE AT START OF NEXT SESSION]
Options:
  A. MILESTONE run — reproduce Wang et al. 2025 HCC key findings using
     GSE166635 through the full Phase 1 pipeline (QC → DEG → GSEA)
  B. scVI batch correction module — deferred from Phase 6, but needed for
     multiome; implement if benchmarking reveals Harmony is insufficient
  C. ADT QC + CLR normalization module (mdata["adt"] path)

Recommended: Option A — close Phase 1 with the benchmark milestone.

ONE goal — do not start scVI, ADT/ATAC modules, or Quarto reports this session.

## Step 1 — Verify all prior tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_deg.py tests/test_gsea.py tests/test_harmony.py \
    tests/test_cluster.py tests/test_pseudobulk_deg.py -v
# Expected: 11 + 8 + 13 + 16 + 14 = 62 passed
```

## Step 2 — Input file for real-data notebook run (Step 10)
```bash
# Input:  data/processed/GSE194122_cite_annotated.h5ad
#         (has obs['cell_type_vote'], obs['batch'], layers['counts'])
# Output: data/processed/GSE194122_cite_pseudobulk_deg.h5ad
# Report: reports/pseudobulk_deg_report.html
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
- pydeseq2 needs min_samples >= 2 per group — filter cell types with < min_samples donors
- pseudobulk test_cell_type_with_too_few_donors_is_skipped uses 3 cell types
  (A+B both pass, C=2 donors is the one that gets skipped) — rest group for A
  must also have >= min_samples pseudo-samples
- pytest.raises(match=...) matches against the exception MESSAGE string,
  not the parameter name — use match="layers\\[" not match="counts_layer"

## Deferred Work (do not start this session)
- scVI batch correction -> deferred to Phase 6 (MultiVI for multiome integration)
- ScType-py + SingleR-py annotation -> docs/ANNOTATION_PLAN.md
- ADT and ATAC QC modules
- Quarto report templates (Phase 2)

## Files Modified Last Session
### New files
- pipeline/modules/downstream/pseudobulk_deg.py   <- new module
- reports/pseudobulk_deg_report.py                <- new report
- tests/test_pseudobulk_deg.py                    <- 14 tests

### Updated files
- notebooks/phase1_qc.ipynb   <- Step 10 added (9 new cells)

### Memory files
- .dev_memory/MODULE_DOCS.md    <- Modules 19-20 added, data flow updated
- .dev_memory/CURRENT_STATUS.md <- Updated
- .dev_memory/NEXT_SESSION.md   <- Updated (this file)
- .dev_memory/PROGRESS.md       <- Updated

## Verify Last Session Works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_pseudobulk_deg.py -v
# Expected: 14 passed
```

## Relevant Context — GSE194122 Data Chain
- data/processed/GSE194122_cite_annotated.h5ad        <- input for pseudobulk_deg()
- data/processed/GSE194122_cite_pseudobulk_deg.h5ad   <- output (new this session)
- data/processed/GSE194122_cite_harmony.h5ad          <- harmony output
- data/processed/GSE194122_cite_harmony_clustered.h5ad <- harmony + clustering output

## Pseudobulk DEG — design notes (for reference)
- pseudobulk_deg() returns pb_dict with keys: results, summary_df, provenance, pairwise, skipped
- pb_dict schema is identical to deg_dict so deg_report.py would work,
  but the dedicated pseudobulk_deg_report.py is preferred (has skipped groups section)
- Provenance stored in uns['omicsage_pseudobulk_deg']
- Report: reports/pseudobulk_deg_report.html
- Notebook Step 10 imports from reports.pseudobulk_deg_report (NOT deg_report)

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, rpy2 3.6.7 (installed but NOT used),
                   pydeseq2 (installed this session — verify version at start)
