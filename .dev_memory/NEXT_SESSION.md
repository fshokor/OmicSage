# NEXT_SESSION.md — Phase 5 Session 1
# Multiome: RNA pipeline run + ATAC QC module
> Updated: June 2026
> Previous phase: Phase 4 (CITE-seq) — COMPLETE

---

## Session Context

**Phase:** 5 — Multiome (RNA + ATAC jointly)
**Session number:** 1 of ~6
**Last thing completed:** Phase 4 closed — cite_corr within-cell-type fix,
five report modules updated, environment.yml and requirements-ci.txt synced,
cite_deg.py n_cells provenance fix. 855+ tests passing.

**File last worked on:** `requirements-ci.txt` / `environment.yml`

---

## Baseline Verification

Run this before writing any new code:

```bash
conda activate omicsage
cd ~/OmicSage
python -m pytest tests/ -q --tb=short 2>&1 | tail -5
```

Expected: 855+ passed, 0 errors. Do not proceed if this fails.

---

## Today's Goals (two tightly scoped deliverables)

### Goal 1 — Run the RNA pipeline on the multiome dataset

The config file `config/runs/GSE194122_multiome.yaml` already exists but
may need updating before running. Verify it, fix any issues, then run
the full RNA pipeline to produce `data/processed/GSE194122_multiome/05_annotated.h5ad`.
This RNA output is the input to all downstream multiome steps.

### Goal 2 — Build `atac_qc.py`

First ATAC-specific module. Takes the raw multiome h5ad, splits out ATAC
peaks, computes per-cell QC metrics, and returns a MuData with `mdata["atac"]`
populated and QC'd.

---

## Step 0 — Before any code: update the config file

The existing `config/runs/GSE194122_multiome.yaml` was written early in Phase 1
for a basic RNA run. It needs to be reviewed and updated to:

1. Point to the correct raw multiome h5ad path
2. Set modality to `multiome` so `qc.py` splits GEX and ATAC correctly
3. Confirm batch key is `batch` (same as CITE — NeurIPS 2021 BMMC uses `batch`)
4. Confirm annotation settings match what Phase 4 used (same donors, same cell types)

**Verify the raw file location first:**
```bash
ls -lh data/benchmark/GSE194122/*multiome* 2>/dev/null || \
ls -lh data/benchmark/GSE194122/ | grep -i multiome
```

**Check what cell types are in the ground truth:**
```python
import scanpy as sc
adata = sc.read_h5ad("data/benchmark/GSE194122/<multiome_file>.h5ad", backed='r')
print(adata.obs['cell_type'].value_counts())
print(adata.var['feature_types'].value_counts())
adata.file.close()
```

Expected output:
- `feature_types`: `Gene Expression` (or `GEX`) + `Peaks` — two modalities in one matrix
- `cell_type`: same 13 NeurIPS cell type labels as the CITE dataset
  (CD4 T, CD8 T, B cell, NK, Monocyte CD14, Monocyte CD16, DC, etc.)

---

## Step 1 — Update `config/runs/GSE194122_multiome.yaml`

Replace the existing file with the template below, adjusting paths as needed.

```yaml
# config/runs/GSE194122_multiome.yaml
# RNA pipeline config for the multiome arm of GSE194122 (NeurIPS 2021 BMMC)
# Run with: python run_pipeline.py --config config/runs/GSE194122_multiome.yaml

dataset_id: GSE194122_multiome
modality: multiome                        # tells qc.py to split GEX + ATAC

input:
  path: data/benchmark/GSE194122/<multiome_filename>.h5ad
  format: h5ad

output:
  processed_dir: data/processed/GSE194122_multiome
  reports_dir: reports/GSE194122_multiome

batch_key: batch
species: human

steps:
  ingest:
    enabled: true

  qc:
    enabled: true
    params:
      min_genes: 200
      max_genes: 7000
      max_mt_pct: 20.0
      min_cells: 3
      run_scrublet: true
      filter_doublets: false            # flag only, same as CITE convention
      exclude_gene_prefixes: []

  normalize:
    enabled: true
    params:
      target_sum: 10000
      n_hvg: 2000
      use_logcounts_layer: true

  reduce:
    enabled: true
    params:
      n_pcs: 50
      n_neighbors: 15
      use_harmony: false                # harmony runs in its own step

  cluster:
    enabled: true
    params:
      resolution: 0.5
      n_iterations: -1

  annotate:
    enabled: true
    params:
      methods: [celltypist, sctype, singler, scanvi]
      vote_weights:
        celltypist: 1
        sctype: 1
        singler: 2
        scanvi: 2
      # Cell type map — same labels as CITE pipeline for cross-modal consistency
      # Fill after inspecting cluster report
      cell_type_map: {}

  harmony:
    enabled: true
    params:
      batch_key: batch
      max_iter_harmony: 20

  deg:
    enabled: true
    params:
      method: wilcoxon
      n_genes: 200
      min_logfc: 0.5
      max_pval_adj: 0.05
      exclude_gene_prefixes: [MT-, RPL, RPS]

  gsea:
    enabled: true
    params:
      gene_sets: [GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022]
      max_pval: 0.05

  pseudobulk:
    enabled: true
    params:
      min_samples: 3
      group_key: cell_type_vote
```

---

## Step 2 — Run the RNA pipeline

```bash
conda activate omicsage
python run_pipeline.py --config config/runs/GSE194122_multiome.yaml \
    2>&1 | tee logs/GSE194122_multiome_rna.log
```

**Expected outputs** in `data/processed/GSE194122_multiome/`:
```
00_ingested.h5ad
01_qc.h5mu          ← MuData (mdata["rna"] + mdata["atac"] both present)
02_normalized.h5ad
03_reduced.h5ad
04_clustered.h5ad
05_annotated.h5ad   ← primary input for multiome steps
06_harmony.h5ad
07_deg.h5ad
08_gsea.h5ad
```

**Key difference from CITE:** `qc.py` already handles multiome — it detects
`feature_types` = GEX + Peaks and returns a MuData with both modalities split.
The RNA pipeline steps (normalize onward) operate on `mdata["rna"]` just as
they do for CITE. The `mdata["atac"]` modality is carried through but untouched
until the multiome phase modules.

**If the pipeline fails at QC with a modality error:**
Check that `modality: multiome` is set in the config. The `qc.py` `_detect_modality()`
function reads `adata.var['feature_types']` — if the column name differs
(e.g. `feature_type` singular, or `ATAC` instead of `Peaks`), print the column
and fix the detection logic before proceeding.

**Disk space warning:** The multiome h5ad is ~2.9 GB. The pipeline produces
~8 checkpoint files. Ensure at least 30 GB free before running:
```bash
df -h ~/OmicSage/data/
```

---

## Step 3 — Build `pipeline/modules/multiome/atac_qc.py`

### Reference to read first (mandatory)
https://www.sc-best-practices.org/chromatin_accessibility/introduction.html
https://www.sc-best-practices.org/chromatin_accessibility/quality_control.html

Key facts to extract from the reference before writing code:
- What are the recommended TSS enrichment score thresholds?
- What does the nucleosome signal distribution look like for good vs bad cells?
- What is the recommended minimum fragment count per cell?
- Does sc-best-practices recommend episcanpy or snapatac2 for peak QC?

### What `atac_qc.py` does

Takes `mdata["atac"]` from the QC output MuData and computes ATAC-specific
per-cell metrics. Returns updated MuData with metrics in `mdata["atac"].obs`.

**ATAC QC metrics to compute:**

| Metric | Description | Typical threshold |
|---|---|---|
| `n_peaks_by_counts` | Number of peaks with > 0 counts per cell | min 500 |
| `total_peak_counts` | Total ATAC counts per cell | min 1000 |
| `pct_counts_in_top_peaks` | % counts in top 500 peaks | < 80% (not too concentrated) |
| `nucleosome_signal` | Ratio of mono-nucleosomal to sub-nucleosomal fragments | < 4 |
| `tss_enrichment` | TSS enrichment score | > 1.0 (loose) / > 2.0 (strict) |
| `atac_predicted_doublet` | From scrublet on peak matrix | flagged only |

**Note on TSS enrichment and nucleosome signal:**
These require the fragment file (`.tsv.gz`), not just the peak count matrix.
For GSE194122 multiome, check whether fragment files are available:
```bash
ls data/benchmark/GSE194122/*fragment* 2>/dev/null
```

If fragment files are NOT available (common for GEO-deposited processed data):
- Skip TSS enrichment and nucleosome signal computation
- Note in provenance: `"fragment_file_available": false`
- Compute only peak-count-based metrics (n_peaks, total_counts, pct_in_top)
- Document this as a known limitation — the multiome report will display a warning

If fragment files ARE available: use `episcanpy` or implement directly via
pysam/tabix. Check sc-best-practices recommendation.

### Function signature

```python
def atac_qc(
    data: MuData,
    min_peaks: int = 500,
    max_peaks: int = 50000,
    min_counts: int = 1000,
    max_nucleosome_signal: float = 4.0,
    min_tss_enrichment: float = 1.0,
    fragment_file: Optional[str] = None,
    run_scrublet: bool = True,
    filter_cells: bool = False,        # flag only by default
    inplace: bool = False,
) -> tuple[MuData, dict]:
```

### Returns `(MuData, atac_qc_dict)`

```python
atac_qc_dict = {
    "metrics": pd.DataFrame,           # per-cell metrics
    "n_cells_before": int,
    "n_cells_after": int,
    "n_peaks": int,
    "fragment_file_available": bool,
    "provenance": dict,                # written to mdata.uns["omicsage_atac_qc"]
}
```

### File location
`pipeline/modules/multiome/atac_qc.py`

Create `pipeline/modules/multiome/__init__.py` as an empty file.

### Tests
`tests/test_atac_qc.py`

Minimum test classes:
- `TestAtacQcInput` — rejects AnnData, requires atac modality in MuData
- `TestAtacQcMetrics` — metrics computed correctly on fixture
- `TestAtacQcFilter` — filter_cells=True removes low-quality cells
- `TestAtacQcProvenance` — uns["omicsage_atac_qc"] written with correct keys
- `TestAtacQcFragmentFallback` — no fragment file → TSS/nucleosome skipped gracefully

Fixture: build a minimal MuData with a synthetic peak count matrix
(n_cells=200, n_peaks=1000) in `mdata["atac"]`. Do NOT use the real 2.9 GB
multiome file in unit tests — it is too slow for CI.

---

## Known issues entering this session

- Stale `__pycache__` from Zone.Identifier files can cause import failures in WSL2.
  If you see unexpected import errors: `find . -name __pycache__ -exec rm -rf {} + 2>/dev/null`
- System pytest at `~/.local/bin/pytest` can shadow the conda env version.
  Always use `python -m pytest`.
- mudata FutureWarnings on `.update()` are harmless — suppress in pipeline entry
  points with `warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")`
- `requirements-ci.txt` and `environment.yml` are now in sync — do not add new
  dependencies to one without updating the other.

---

## Files to create this session

```
pipeline/modules/multiome/__init__.py       (empty)
pipeline/modules/multiome/atac_qc.py        (new module)
tests/test_atac_qc.py                       (new tests)
config/runs/GSE194122_multiome.yaml         (update existing)
logs/GSE194122_multiome_rna.log             (generated by pipeline run)
```

## Files that may need updating

```
pipeline/modules/qc/qc.py     — only if modality detection fails on multiome file
reports/GSE194122_multiome/   — generated by RNA pipeline run
```

---

## End-of-session checklist

- [ ] RNA pipeline ran to completion for GSE194122_multiome
- [ ] `data/processed/GSE194122_multiome/05_annotated.h5ad` exists
- [ ] `atac_qc.py` written and returns `(MuData, dict)`
- [ ] `uns["omicsage_atac_qc"]` written with provenance
- [ ] All new tests pass: `python -m pytest tests/test_atac_qc.py -v`
- [ ] Full test suite still passes: `python -m pytest tests/ -q | tail -5`
- [ ] `NEXT_SESSION.md` updated for Session 2 (ATAC reduce + LSI)
- [ ] `PROGRESS.md` updated: Phase 5 Session 1 ticked

---

## Phase 5 roadmap (all sessions)

| Session | Primary deliverable | Secondary |
|---|---|---|
| **1 (this)** | RNA pipeline run on multiome data + `atac_qc.py` | Config file update |
| 2 | `atac_reduce.py` — TF-IDF → LSI → UMAP | `atac_qc_report.py` |
| 3 | `atac_annotate.py` — gene activity scores + cluster annotation | `atac_reduce_report.py` |
| 4 | `multiome_integration.py` — MultiVI (scvi-tools) joint embedding | `atac_annotate_report.py` |
| 5 | `multiome_deg.py` — joint RNA DEG + differential accessibility peaks | `multiome_integration_report.py` |
| 6 | `run_multiome_pipeline.py` + combined report | All remaining report modules |

---

## Reference URLs for this session

- sc-best-practices ATAC QC: https://www.sc-best-practices.org/chromatin_accessibility/quality_control.html
- sc-best-practices ATAC intro: https://www.sc-best-practices.org/chromatin_accessibility/introduction.html
- episcanpy docs: https://episcanpy.readthedocs.io
- snapatac2 docs: https://kzhang.org/SnapATAC2
- MultiVI paper (scvi-tools): https://www.nature.com/articles/s41592-023-01945-9
- NeurIPS 2021 BMMC dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122

---

## Key architectural decisions for Phase 5

**Why MultiVI for integration (not WNN)?**
WNN via `muon.pp.neighbors` hangs on small fixtures (known pynndescent issue,
flagged in Phase 4 Session 6). MultiVI is pure Python, handles RNA + ATAC
jointly in a single VAE, and is the sc-best-practices recommended method for
multiome integration. WNN will be revisited in Phase 6 (Multiome extended) if
needed.

**Namespace conventions (do not change):**
- `mdata["rna"]` — RNA modality throughout (matches CITE convention)
- `mdata["atac"]` — ATAC modality
- `mdata["rna"].obsm["X_pca"]` — RNA PCA
- `mdata["atac"].obsm["X_lsi"]` — ATAC LSI (not PCA — LSI is the correct
  dimensionality reduction for sparse binary peak matrices)
- `mdata.obsm["X_multivi"]` — joint MultiVI embedding
- `mdata["atac"].obs["atac_celltype"]` — ATAC-derived annotation
- `mdata["rna"].obs["cell_type_vote"]` — RNA consensus annotation (from Phase 1)
- Batch key: `batch` throughout (same as CITE, same donors)
