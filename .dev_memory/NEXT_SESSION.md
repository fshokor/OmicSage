# NEXT_SESSION.md — Phase 7 Session 5
# Spatial Transcriptomics: Downstream Analysis
> Updated: June 2026
> Previous session: Phase 7 Sessions 2–4 COMPLETE (reduce, cluster, deconvolve, reports, runner)

---

## Session Context

**Phase:** 7 — Spatial Transcriptomics
**Session number:** 5 of 5 (final Phase 7 session)
**Last things completed (Sessions 2–4):**
  - `spatial_reduce.py`: normalize, HVG, PCA, spatial neighbours graph
  - `spatial_cluster.py`: Leiden clustering + Moran's I SVGs
  - `spatial_deconvolve.py`: cell2location deconvolution (graceful skip when no ref)
  - `spatial_ingest.py` updated: ENSEMBL ID swap, layers["counts"] preservation, MT gene stripping
  - All 4 spatial report files rewritten to match RNA `_render_page` / `<main>/<section>` pattern
  - `spatial_combined_report.py`: tabbed combined report (identical UI to RNA combined_report.py)
  - `run_spatial_pipeline.py`: full runner with --step / --from-step / --to-step / --force / checkpointing
  - `config/runs/kuppe_heart.yaml`: Kuppe 2022 human heart config
  - `tests/test_spatial_reduce.py`, `tests/test_spatial_cluster.py`, `tests/test_spatial_deconvolve.py`

**Benchmark dataset:** Kuppe et al. 2022 human myocardial infarction (control samples)
  - Visium: `data/benchmark/kuppe_visium_human_heart_2022_control.h5ad`
  - snRNA-seq reference: `data/benchmark/kuppe_snRNA_human_heart_2022_control.h5ad`
  - cell_type_key: `cell_type_original`, batch_key_ref: `donor_id`, batch_key_st: `patient`

**Current test count:** 1378 passed, 57 skipped (verified end of Session 4)

---

## Baseline Verification

```bash
conda activate omicsage
cd ~/OmicSage
python -m pytest tests/ -q --tb=short 2>&1 | tail -5
```

Expected: 1378 passed, 57 skipped. Confirm before writing any new code.

Also verify spatial imports:
```bash
python -c "
from pipeline.modules.spatial.spatial_ingest import spatial_ingest
from pipeline.modules.spatial.spatial_reduce import spatial_reduce
from pipeline.modules.spatial.spatial_cluster import spatial_cluster
from pipeline.modules.spatial.spatial_deconvolve import spatial_deconvolve
print('all spatial imports OK')
"
```

---

## Today's Goal — `spatial_downstream.py`

**Single deliverable:** one module covering all spatial downstream analyses,
plus its report (`spatial_downstream_report.py`), tests, and runner update.

This is the final Phase 7 session. After this, Phase 7 is complete.

---

## Analyses to implement (all in one file)

### 1. Cell-type specific gene expression in spatial context
- For each cell type from deconvolution, identify genes whose expression
  correlates with that cell type's abundance across spots
- Method: Spearman correlation between `obs[cell_type]` abundance and
  gene expression across spots
- Output: `uns["celltype_marker_genes"]` — dict of cell_type → top 20 correlated genes
- Requires: `obsm["q05_cell_abundance_w_sf"]` (from deconvolution)
- Graceful skip if deconvolution was skipped

### 2. Spatially variable genes per cell type (cell-type-specific SVGs)
- Filter Moran's I results to spots enriched for each cell type
  (spots where cell type abundance > median abundance for that cell type)
- Run `sq.gr.spatial_autocorr(mode="moran")` on the enriched subset
- Output: `uns["celltype_svg"]` — dict of cell_type → DataFrame of SVGs
- Requires: `uns["moranI"]` (from clustering) + deconvolution abundances
- Graceful skip if either is absent

### 3. Spatial co-occurrence (`sq.gr.co_occurrence`)
- Which cell types tend to co-localise spatially at different distance ranges
- Output: stored in `uns["co_occurrence"]` by squidpy
- Report: heatmap of co-occurrence scores at a reference distance
- Requires: deconvolution cell type labels in `obs["dominant_cell_type"]`
- Graceful skip if dominant_cell_type absent

### 4. Cell type neighbourhood enrichment (`sq.gr.nhood_enrichment`)
- Which cell types are enriched in each other's neighbourhoods
- Method: permutation test on spatial adjacency graph
- Output: stored in `uns["nhood_enrichment"]` by squidpy
- Report: heatmap of enrichment z-scores
- Requires: `obs["dominant_cell_type"]` + `obsp["spatial_connectivities"]`

### 5. Ligand-receptor communication (`sq.gr.ligrec`)
- Which cell type pairs are likely communicating based on proximity + expression
- Uses CellChatDB (human) via squidpy's built-in database
- Output: stored in `uns["ligrec"]` by squidpy
- Report: dot plot of top LR pairs by cell type pair
- Requires: deconvolution cell type labels
- Graceful skip if cell types absent

### 6. Pathway enrichment on top SVGs (reuse gseapy)
- Run GSEA on the ranked SVG list from Moran's I
- Gene sets: GO Biological Process (human) via gseapy
- Output: `uns["svg_gsea"]` — top enriched pathways
- Requires: `uns["moranI"]`
- Graceful skip if moranI absent

---

## Function signature

```python
def spatial_downstream(
    adata: AnnData,
    # Cell-type specific expression
    run_celltype_expression: bool = True,
    n_marker_genes: int = 20,
    # Cell-type specific SVGs
    run_celltype_svg: bool = True,
    svg_n_genes: Optional[int] = None,
    # Co-occurrence
    run_co_occurrence: bool = True,
    co_occurrence_interval: Optional[list] = None,  # default: squidpy auto
    # Neighbourhood enrichment
    run_nhood_enrichment: bool = True,
    n_perms_nhood: int = 1000,
    # Ligand-receptor
    run_ligrec: bool = True,
    ligrec_n_perms: int = 1000,
    ligrec_organism: str = "human",      # "human" or "mouse"
    # Pathway enrichment on SVGs
    run_svg_gsea: bool = True,
    svg_gsea_gene_sets: str = "GO_Biological_Process_2023",
    svg_gsea_organism: str = "Human",
    # Cell type column (from deconvolution)
    dominant_celltype_key: str = "dominant_cell_type",
    # General
    n_jobs: int = 1,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
```

All `run_*` flags are True by default but each analysis gracefully skips
if its required inputs are absent (moranI, deconvolution abundances, etc.)
and records `skipped=True` in provenance.

---

## Output keys written to AnnData

```
uns["celltype_marker_genes"]     dict: cell_type → list of top correlated genes
uns["celltype_svg"]              dict: cell_type → DataFrame (Moran's I per subset)
uns["co_occurrence"]             set by sq.gr.co_occurrence
uns["nhood_enrichment"]          set by sq.gr.nhood_enrichment
uns["ligrec"]                    set by sq.gr.ligrec
uns["svg_gsea"]                  DataFrame: top enriched pathways from SVG list
uns["omicsage_spatial_downstream"]  provenance dict
```

---

## Report content (`spatial_downstream_report.py`)

Sections using the standard `_render_page` / `<main>/<section>` pattern:

1. **Run Summary** — stat cards: n analyses run, n cell types, n SVGs tested
2. **Cell-type marker genes** — table: top 10 correlated genes per cell type
3. **Cell-type specific SVGs** — table: top 5 SVGs per cell type
4. **Co-occurrence** — heatmap: cell type co-occurrence at reference distance
5. **Neighbourhood enrichment** — heatmap: enrichment z-scores between cell types
6. **Ligand-receptor interactions** — dot plot: top LR pairs
7. **SVG pathway enrichment** — bar chart: top 10 enriched GO terms

All sections gracefully show "not run / data not available" when skipped.

---

## Files to create this session

```
pipeline/modules/spatial/spatial_downstream.py          ← new
reports/templates/spatial/spatial_downstream_report.py  ← new
tests/test_spatial_downstream.py                        ← new
```

Update:
```
run_spatial_pipeline.py   ← add step 6: downstream
reports/spatial_combined_report.py  ← add tab for downstream report
config/runs/kuppe_heart.yaml        ← add downstream section
```

---

## Pre-implementation research required

Before writing any code, fetch and read:

1. **sq.gr.co_occurrence API**
   https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.co_occurrence.html

2. **sq.gr.nhood_enrichment API**
   https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.nhood_enrichment.html

3. **sq.gr.ligrec API**
   https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.ligrec.html

4. **sc-best-practices spatial downstream chapter**
   https://www.sc-best-practices.org/spatial/spatially_variable_genes.html

---

## Runner update

Add step 6 to `run_spatial_pipeline.py`:

```
STEP_ORDER = ["ingest", "qc", "reduce", "cluster", "deconvolve", "downstream"]
STEP_OUTPUT["downstream"] = "06_downstream.h5ad"
STEP_REPORT["downstream"] = "spatial_downstream_report.html"
STEP_PREDECESSOR["downstream"] = "deconvolve"
```

Config section to add to `kuppe_heart.yaml`:
```yaml
  downstream:
    run_celltype_expression: true
    n_marker_genes: 20
    run_celltype_svg: true
    run_co_occurrence: true
    run_nhood_enrichment: true
    run_ligrec: true
    ligrec_organism: "human"
    run_svg_gsea: true
    svg_gsea_gene_sets: "GO_Biological_Process_2023"
    svg_gsea_organism: "Human"
    n_jobs: 4
```

---

## Combined report update

Add to `TAB_REGISTRY` in `spatial_combined_report.py`:
```python
("spatial_downstream_report.html", "Downstream", "🔗"),
```

---

## Key implementation notes

- **sq.gr.co_occurrence** requires a cluster key (not continuous abundances).
  Use `obs["dominant_cell_type"]` — the argmax cell type per spot written
  by `spatial_deconvolve_report.py`. If absent, skip gracefully.

- **sq.gr.nhood_enrichment** also requires a cluster key — same key.

- **sq.gr.ligrec** requires gene symbols in var_names but after ingest
  we swap var_names to ENSEMBL IDs. Use `var["feature_name"]` as the
  gene symbol source and temporarily rename before ligrec call,
  or check if `var["feature_name"]` exists and swap back.

- **Cell-type-specific SVGs**: create a boolean mask per cell type
  (`adata.obs[ct] > adata.obs[ct].median()`), subset adata to those spots,
  run `sq.gr.spatial_autocorr` on the subset. The spatial neighbours graph
  must be recomputed for the subset (smaller graph). Use n_perms=None
  (analytical p-values) for speed.

- **GSEA on SVGs**: the SVG gene list from `uns["moranI"]` is already
  ranked by Moran's I score. Use `gseapy.prerank` with the ranked gene
  list. Gene names must be symbols — same ENSEMBL/symbol swap issue as ligrec.

- **Never densify** large sparse matrices. Cell-type expression correlation
  uses column slicing on CSC format.

---

## DECISIONS.md entries to add this session

```markdown
### Downstream: dominant_cell_type key
- sq.gr.co_occurrence and sq.gr.nhood_enrichment require a categorical
  cluster key, not continuous abundances.
- We use obs["dominant_cell_type"] = argmax of q05_cell_abundance_w_sf.
- This is written during spatial_deconvolve_report figure generation
  (side effect). We promote it to a first-class output of spatial_deconvolve.py.

### Downstream: gene symbol vs ENSEMBL ID
- After ingest, var_names are ENSEMBL IDs; var["feature_name"] has symbols.
- ligrec and GSEA require gene symbols.
- Solution: pass genes=adata.var["feature_name"].tolist() to ligrec,
  and build the GSEA ranked dict using feature_name as index.

### Downstream: cell-type SVGs use analytical Moran's I
- n_perms=None (analytical p-values) for speed.
- Permutation testing available but slow; documented in config.
```

---

## Phase 7 completion checklist

After this session, Phase 7 is complete. Tick these off before moving to Phase 8:

- [ ] All 6 spatial steps run end-to-end on Kuppe data
- [ ] Combined report shows all 5 tabs (QC, Reduce, Cluster, Deconvolve, Downstream)
- [ ] cell2location full run completed (250 epochs ref + 30k epochs spatial)
- [ ] Test count: expected ~1378 + ~45 new = ~1423 passed
- [ ] `PROGRESS.md` updated: Phase 7 complete
- [ ] `CURRENT_STATUS.md` updated
- [ ] `DECISIONS.md` updated with all Phase 7 decisions

---

## Phase 8 preview (Streamlit no-code UI)

After Phase 7, the next phase is the Streamlit interface. Key deliverables:
- File upload (h5ad, Visium directory)
- Config builder (guided form)
- Pipeline progress tracking
- Report viewer (embedded HTML)
- Reference: `ui/app.py` stub already exists in repo
