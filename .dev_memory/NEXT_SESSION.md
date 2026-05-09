## Session Context
Date: 2026-05-09
Phase: 1 — Core scRNA Pipeline
Last thing completed: DEG module — deg.py, test_deg.py (11 tests), deg_report.py,
                      notebook Step 6, MODULE_DOCS.md updated
File last worked on: notebooks/phase1_qc.ipynb (Step 6 DEG section added)

## Today's Goal
GSEA module — implement gsea.py (GO/KEGG/Reactome pathway enrichment
from deg_dict['results'], HTML report, notebook Step 7).
ONE goal — do not start batch correction or annotation improvements this session.

## Step 1 — Verify DEG module still works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_deg.py -v
# Expected: 11 passed
```

## Step 2 — Implement gsea.py
File to create: pipeline/modules/qc/gsea.py

GSEA requirements:
- Input: deg_dict['results'] (dict of {group: DataFrame}) + adata (for gene universe)
- Method: over-representation analysis (ORA) via gseapy (Fisher exact + BH correction)
- Gene sets: GO Biological Process, KEGG, Reactome (configurable)
- Per-group: top N significant DEGs (up-regulated only by default) as query gene list
- Gene universe: all genes in adata.var_names
- Output: gsea_dict with {group: DataFrame(Term, Overlap, P-value, Adjusted P-value, Genes)}
- Write provenance to adata.uns['omicsage_gsea']
- Graceful handling: skip groups with < min_genes significant DEGs (warn, don't crash)

Proposed API:
```python
from pipeline.modules.qc.gsea import gsea

adata_gsea, gsea_dict = gsea(
    adata_deg,
    deg_dict=deg_dict,
    gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"],
    min_logfc=0.25,          # re-filter DEGs before enrichment (or pass pre-filtered)
    max_pval_adj=0.05,
    top_n_genes=None,        # None = use all significant DEGs per group
    min_genes=5,             # skip group if fewer than 5 DEGs pass filters
    organism="human",        # "human" | "mouse"
    inplace=False,
)
```

gsea_dict keys:
  'results'    — {group: DataFrame(Term, Overlap, P-value, Adjusted P-value, Genes)}
  'summary_df' — top 3 terms per group, long format
  'provenance' — same as uns['omicsage_gsea']

uns['omicsage_gsea'] keys:
  groupby, gene_sets, organism, min_logfc, max_pval_adj,
  top_n_genes, n_groups_tested, n_groups_skipped,
  gseapy_version, omicsage_module, timestamp

## Step 3 — Write test_gsea.py
File to create: tests/test_gsea.py

Tests to write:
- test_gsea_returns_anndata_and_dict()
- test_gsea_uns_provenance_keys()
- test_gsea_output_columns()          — Term, Overlap, P-value, Adjusted P-value, Genes
- test_gsea_pval_range()              — all pvals in [0, 1]
- test_gsea_every_group_has_results() — or is in skipped list
- test_gsea_min_genes_skip()          — groups with <min_genes DEGs skipped with warning
- test_gsea_inplace_false()           — caller object unchanged
- test_gsea_gene_sets_param()         — single gene set string also accepted

## Step 4 — Implement gsea_report.py
File to create: reports/gsea_report.py

Report sections:
- Summary: groups tested, gene sets queried, total significant pathways
- Top pathways per group table (top 5, adj. p-value, overlap, gene list)
- Bar chart per group: top 10 pathways by -log10(adj. p-value)
- Bubble plot across groups (optional, if ≤9 groups): pathway × group, size = overlap, colour = adj.p
- HTML output matching deg_report.py style

## Step 5 — Update notebook
File to edit: notebooks/phase1_qc.ipynb

Add Step 7 section:
- Load data/processed/GSE194122_cite_deg.h5ad (and deg_dict, or re-run deg())
- Run gsea() — GO BP + KEGG + Reactome
- Sanity check: known immune pathways appear for expected cell types
  (e.g. T cell activation GO terms in T cells, phagocytosis in Monocytes)
- Generate HTML report
- Save data/processed/GSE194122_cite_gsea.h5ad

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- OldFormatWarning from GSE194122 — harmless
- Notebook must be opened from OmicSage root
- seurat_v3 HVG flavor needs ≥500 cells
- numpy bool: use bool() cast before `is True` checks in tests
- float dict keys must be stringified before storing in uns
- obs['cell_type_vote'] is the consensus column (not obs['cell_type'])
- obs['cell_type_groundtruth'] holds the preserved publication ground-truth
- rpy2 NOT used — do not attempt R integration
- gseapy dependency: add to requirements-ci.txt and environment.yml if not present

## Deferred Work (do not start this session)
- ScType-py + SingleR-py annotation → docs/ANNOTATION_PLAN.md
- Harmony + scVI batch correction
- ADT and ATAC QC modules
- Pseudobulk DEG (DESeq2-style)

## Files Modified Last Session
- pipeline/modules/qc/deg.py              ← CREATED
- tests/test_deg.py                       ← CREATED
- reports/deg_report.py                   ← CREATED
- notebooks/phase1_qc.ipynb              ← Step 6 DEG section added (9 cells)
- .dev_memory/MODULE_DOCS.md             ← Modules 11–14 added
- .dev_memory/CURRENT_STATUS.md          ← Updated
- .dev_memory/NEXT_SESSION.md            ← Updated (this file)
- .dev_memory/PROGRESS.md               ← Updated

## Verify Last Session Works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_deg.py -v
# Expected: 11 passed
```

## Relevant Context — GSE194122 Data Chain
- data/processed/GSE194122_cite_annotated.h5ad   ← input for deg()
- data/processed/GSE194122_cite_deg.h5ad         ← input for gsea()
- data/processed/GSE194122_cite_gsea.h5ad        ← output of next session

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   rpy2 3.6.7 (installed but NOT used)
Check gseapy before starting: pip show gseapy
If missing: pip install gseapy --break-system-packages (or conda install -c bioconda gseapy)
