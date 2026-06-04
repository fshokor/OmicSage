# NEXT_SESSION.md — Phase 5 Session 3
# Multiome: atac_annotate.py — gene activity scores + label transfer
> Updated: June 2026
> Previous session: Phase 5 Session 2 — COMPLETE

---

## Session Context

**Phase:** 5 — Multiome (RNA + ATAC jointly)
**Session number:** 3 of ~6
**Last thing completed:** Phase 5 Session 2 closed — `atac_reduce.py` (TF-IDF → LSI →
UMAP → Leiden), `atac_reduce_report.py`, `tests/test_atac_reduce.py`. All tests passing.

**Files created last session:**
```
pipeline/modules/multiome/atac_reduce.py
reports/templates/multiome/atac_reduce_report.py
tests/test_atac_reduce.py
```

**Files created Session 1 (for reference):**
```
pipeline/modules/multiome/__init__.py       (empty)
pipeline/modules/multiome/atac_qc.py
reports/templates/multiome/atac_qc_report.py
tests/test_atac_qc.py
run_pipeline.py                             (ATAC save fix — 01_qc_atac.h5ad)
config/runs/GSE194122_multiome.yaml         (synced to CITE config)
```

---

## Baseline Verification

Run this before writing any new code:

```bash
conda activate omicsage
cd ~/OmicSage
python -m pytest tests/ -q --tb=short 2>&1 | tail -5
```

Expected: 855+ passed (+ new atac_qc and atac_reduce tests), 0 errors.
Do not proceed if this fails.

---

## Today's Goal — `atac_annotate.py`

Build the ATAC annotation module. Two steps:

### Step A — Gene activity scores
Compute a proxy gene expression matrix from peak counts by summing peaks
that overlap each gene's body + 2kb upstream promoter window.

Source for peak-to-gene mapping: `mdata["atac"].uns["atac"]["peak_annotation"]`
This is a DataFrame with columns `peak`, `gene_name`, `distance`, `peak_type`
populated automatically by muon when loading 10x h5. It already exists in the
`01_qc_atac.h5ad` checkpoint.

Store result in: `mdata["atac"].obsm["gene_activity"]`
(shape: n_cells × n_genes; dense numpy array or sparse matrix)
Also store gene names in: `mdata["atac"].uns["gene_activity_var_names"]`

**Important:** if `peak_annotation` is absent (e.g. in unit test fixtures),
fall back to a proximity-based approach using peak var_names parsed as
`chr:start-end` coordinates. Document this clearly.

### Step B — Label transfer from RNA
Transfer `cell_type_vote` labels from the RNA AnnData to ATAC cells using
majority vote per Leiden cluster (`atac_leiden`).

Input: `atac.obs["atac_leiden"]` + RNA AnnData with `obs["cell_type_vote"]`
Method: for each `atac_leiden` cluster, find the majority `cell_type_vote`
among cells in that cluster (barcodes must match between RNA and ATAC).
Output: `atac.obs["atac_celltype"]`

This is simpler and more robust than KNN label transfer for multiome data
because we already have a 1:1 cell correspondence (same barcodes).
Reference: Seurat v5 and Signac vignettes both use label transfer for
unpaired data; for paired multiome data, direct barcode matching is preferred.

---

## Function signature

```python
def atac_annotate(
    atac: AnnData,
    rna: Optional[AnnData] = None,        # RNA AnnData with obs["cell_type_vote"]
    peak_annotation: Optional[pd.DataFrame] = None,  # override uns["atac"]["peak_annotation"]
    promoter_upstream_bp: int = 2000,     # upstream window for gene activity
    min_peaks_per_gene: int = 1,          # min peaks overlapping a gene to include it
    leiden_key: str = "atac_leiden",      # cluster key for label transfer
    rna_label_key: str = "cell_type_vote",# RNA obs column to transfer
    inplace: bool = False,
) -> tuple[AnnData, dict]:
```

Returns `(AnnData, dict)` with:
- `atac.obsm["gene_activity"]`              gene activity score matrix
- `atac.uns["gene_activity_var_names"]`     gene names (list)
- `atac.obs["atac_celltype"]`              transferred label (or "Unknown" if no RNA)
- `atac.uns["omicsage_atac_annotate"]`     provenance

---

## Files to create this session

```
pipeline/modules/multiome/atac_annotate.py        (new module)
reports/templates/multiome/atac_annotate_report.py (new report)
tests/test_atac_annotate.py                        (new tests)
```

---

## Files to upload at session start

Claude needs these two files to match conventions:
- `pipeline/modules/cite/adt_annotate.py`   — to match (AnnData, dict) return pattern
- `tests/test_adt_annotate.py`              — to match fixture and assertion style

---

## ATAC data reminder

The actual GSE194122 multiome ATAC AnnData has:
```
AnnData object with n_obs × n_vars = 68715 × 116490
obs: 'cell_type', 'GEX_pct_counts_mt', 'batch', 'DonorID', 'Site',
     'Samplename', 'DonorNumber', 'DonorAge', 'DonorBMI', 'DonorBloodType',
     'DonorRace', 'Ethnicity', 'DonorGender', 'DonorSmoker', 'Modality',
     'ATAC_nCount_peaks', 'ATAC_atac_fragments', 'ATAC_reads_in_peaks_frac',
     'ATAC_blacklist_fraction', 'ATAC_nucleosome_signal', 'sample'
var: 'feature_types', 'gene_id'
uns: 'omicsage_source'
```

- `cell_type` in obs → preserved as `cell_type_groundtruth` by `atac_qc.py`
- Fragment file NOT available — CellRanger-ARC pre-computed metrics used
- Batch key: `batch` throughout (same as CITE)
- Peak annotation may need to be loaded separately from the original h5 file
  if `uns["atac"]["peak_annotation"]` is absent after round-tripping through h5ad

---

## Key architectural decisions for this session

**Why majority vote per cluster (not KNN)?**
For multiome data, RNA and ATAC are measured in the same cells — barcodes match
exactly. The simplest and most reliable approach is: for each ATAC Leiden cluster,
find the majority `cell_type_vote` among the RNA cells with matching barcodes.
KNN label transfer (Seurat/Signac approach) is needed for unpaired data only.

**Namespace conventions (do not change):**
- `atac.obsm["gene_activity"]` — gene activity matrix (n_cells × n_genes)
- `atac.uns["gene_activity_var_names"]` — gene names corresponding to gene_activity columns
- `atac.obs["atac_celltype"]` — transferred label (namespaced, no collision with RNA)
- `atac.obs["cell_type_groundtruth"]` — NeurIPS ground truth (set by atac_qc.py)
- `atac.obs["atac_leiden"]` — ATAC Leiden clusters (set by atac_reduce.py)

---

## Known issues / environment quirks

- Stale `__pycache__` from Zone.Identifier files can cause import failures in WSL2.
  If unexpected import errors: `find . -name __pycache__ -exec rm -rf {} + 2>/dev/null`
- Always use `python -m pytest` (not `pytest`) to avoid conda env shadowing.
- mudata FutureWarnings on `.update()` are harmless.
- `requirements-ci.txt` and `environment.yml` must be kept in sync.

---

## Phase 5 roadmap (all sessions)

| Session | Primary deliverable | Status |
|---|---|---|
| 1 | `atac_qc.py` + report + tests + pipeline fix | ✅ DONE |
| 2 | `atac_reduce.py` (TF-IDF → LSI → UMAP) + report + tests | ✅ DONE |
| **3 (this)** | **`atac_annotate.py` — gene activity + label transfer** | ⬜ |
| 4 | `multiome_integration.py` — MultiVI joint embedding | ⬜ |
| 5 | `multiome_deg.py` — RNA DEG + differential accessibility | ⬜ |
| 6 | `run_multiome_pipeline.py` + combined report | ⬜ |

---

## Reference URLs for this session

- sc-best-practices ATAC QC: https://www.sc-best-practices.org/chromatin_accessibility/quality_control.html
- sc-best-practices paired integration (MultiVI section): https://www.sc-best-practices.org/multimodal_integration/paired_integration.html#multiome-data
- Signac PBMC vignette (GeneActivity function): https://stuartlab.org/signac/articles/pbmc_vignette.html
- Seurat v5 ATAC integration (GeneActivity + label transfer): https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
- NeurIPS 2021 BMMC dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122
