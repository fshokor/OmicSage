# OmicSage — Current Status
> Last updated: June 2026

## Phase
Phase 5 — Multiome (RNA + ATAC jointly)
Session 1 planned. Phase 4 (CITE-seq) fully complete.

## Test Count
855+ passing, 1 skipped (as of end of Phase 4)

---

## ✅ Phase 0 — Foundation (COMPLETE)
- Repo structure scaffolded
- Docker base images (Python + R) — defined, not yet built locally
- CI/CD via GitHub Actions — configured
- YAML config schema — defined
- .dev_memory/ system — initialized

---

## ✅ Phase 1 — Core scRNA Pipeline (COMPLETE)

### Ingestion (pipeline/modules/qc/ingest.py)
- Auto-detects 10x MEX, H5, AnnData formats
- Moves normalized values out of .X, puts raw counts in .X
- Handles GSE194122 CITE-seq and multiome, GSE166635 HCC
- load_dataset_dir() — scans parent folder, loads all MTX subfolders,
  concatenates with obs['sample'] + obs['batch'] per subfolder name

### QC (pipeline/modules/qc/qc.py)
- Modality-aware: auto-detects GEX / ADT / ATAC from var['feature_types']
- Returns MuData with mdata["rna"], mdata["adt"], mdata["atac"] slots
- MT%, ribo%, hb%, genes/cell, Scrublet doublets
- 56 tests passing in tests/test_qc.py

### Normalization (pipeline/modules/qc/normalize.py)
- CP10K → log1p, layers['counts'] + layers['logcounts']
- HVG selection — top 2000, seurat_v3 flavor, batch_key support
- 12 tests passing in tests/test_normalize.py

### Dimensionality Reduction (pipeline/modules/qc/reduce.py)
- PCA (HVG subset, n_comps=50, arpack), elbow detection via kneed
- Neighbor graph, UMAP (always), t-SNE (optional)
- 12 tests passing in tests/test_reduce.py

### Leiden Clustering (pipeline/modules/clustering/cluster.py)
- Resolution sweep, silhouette scoring, best_resolution_override
- obs['leiden'] + obs['leiden_harmony'] convenience keys
- compute_ari() for cross-annotation comparison
- 16 tests passing in tests/test_cluster.py

### Cell-Type Annotation (pipeline/modules/annotation/annotate.py)
- CellTypist + ScType + SingleR + scANVI — 4-way weighted vote
- obs['cell_type_vote'], obs['cell_type_confidence']
- 18 tests passing, 1 skipped in tests/test_annotate.py

### DEG (pipeline/modules/downstream/deg.py)
- Wilcoxon rank-sum, one-vs-rest, BH correction
- exclude_gene_prefixes param
- 11 tests passing in tests/test_deg.py

### GSEA (pipeline/modules/downstream/gsea.py)
- ORA via gseapy.enrichr — GO BP / KEGG / Reactome
- direction param: up / down / both
- 8 tests passing in tests/test_gsea.py (Enrichr mocked)

### Harmony Batch Correction (pipeline/modules/integration/harmony_correct.py)
- obsm['X_pca_harmony'], obsm['X_umap_harmony']
- neighbors_harmony graph in uns + obsp
- 13 tests passing in tests/test_harmony.py

### Pseudobulk DEG (pipeline/modules/downstream/pseudobulk_deg.py)
- DESeq2 Wald tests via pydeseq2, one-vs-rest per cell type
- min_cells + min_samples filters with graceful skip
- 14 tests passing in tests/test_pseudobulk_deg.py

### Reports (reports/)
- normalization_report.py, reduce_report.py, cluster_report.py,
  annotate_report.py, deg_report.py, gsea_report.py,
  harmony_report.py, pseudobulk_deg_report.py
- combined_report.py → 00_combined_report.html (tabbed)

### Pipeline Runner (run_pipeline.py)
- Config-driven, --from-step / --to-step / --step / --force flags
- Step-level .h5ad caching, startup validation

### Config System (config/)
- config/schema.yaml — platform master schema
- config/runs/GSE194122.yaml
- config/runs/GSE166635.yaml
- config/runs/GSE194122_multiome.yaml ← needs update for Phase 5

### MILESTONE: Wang et al. 2025 HCC Benchmark ✓
- Full pipeline on GSE166635, key markers and pathways recovered

---

## ✅ Phase 2 — Report Engine (COMPLETE)
- Per-step HTML reports for all Phase 1 steps
- combined_report.py producing tabbed 00_combined_report.html
- MILESTONE: biologist receives full report from one command ✓

---

## ✅ Phase 3 — AI Layer (COMPLETE / ⏸ PAUSED BY DECISION)
Decision (2026-05-17): AI layer fully built and tested. ai_features: false
is the default going forward. Code is NOT deleted. See DECISIONS.md.

Modules built (all tested):
- ai/pipeline_advisor.py (A1) — 13 tests
- ai/clustering_advisor.py (A2) — 22 tests
- ai/cluster_annotator.py (A3) — 23 tests
- ai/deg_validator.py (B1) — 25 tests
- ai/coherence_reviewer.py (B2) — 22 tests
- ai/downstream_suggester.py (B3) — 20 tests
- ai/narrative_generator.py (C1) — 18 tests
- ai/report_writer.py (C2) — 20 tests
- ai/report_reviewer.py (D1) — 18 tests
- Infrastructure (_base, _config_gate, _audit_log, _llm_client,
  _skill_loader) — 20 tests

Manual review mode: CITE_MASTER_PROMPT.md written (June 2026).
Attach + send workflow — no infrastructure required.

---

## ✅ Phase 4 — CITE-seq Pipeline (COMPLETE)

### ADT Normalize (pipeline/modules/cite/adt_normalize.py)
- CLR normalization via muon (clr_axis=0, per-protein)
- DSB deferred (no empty droplets in GSE194122)
- layers['adt_clr'] written

### ADT Doublets (pipeline/modules/cite/adt_doublets.py)
- Marker-pair scoring (CD3/CD19, CD3/CD14)
- obs['adt_predicted_doublet'] — flagged only by default

### ADT Reduce (pipeline/modules/cite/adt_reduce.py)
- PCA + neighbors + UMAP on CLR matrix
- obsm keys: X_pca_adt, X_umap_adt

### ADT Harmony (pipeline/modules/cite/adt_harmony.py)
- harmonypy called directly (not scanpy wrapper)
- Runtime shape detection: if z.shape[0] != n_obs: z = z.T
- obsm['X_pca_harmony_adt']

### ADT Annotate (pipeline/modules/cite/adt_annotate.py)
- Leiden clustering on ADT embedding
- obs['adt_celltype'], obs['adt_celltype_manual']
- Uses adt_celltype key (not cell_type_vote) to avoid collision

### CITE Integration (pipeline/modules/cite/cite_integration.py)
- MOFA+ (via muon) + totalVI (via scvi-tools)
- WNN deferred — pynndescent hang on small fixtures (@pytest.mark.slow)
- scib metrics: iLISI, graph_conn, cLISI, ASW_label for both methods
- obsm['X_mofa'], obsm['X_totalvi']

### CITE DEG / DPE (pipeline/modules/cite/cite_deg.py)
- Wilcoxon DPE on CLR ADT values (one-vs-rest)
- Cross-modal RNA DEG using ADT-defined labels
- n_cells now written to provenance (June 2026 fix)

### CITE GSEA (pipeline/modules/cite/cite_gsea.py)
- ORA on RNA arm using ADT-defined cell type labels

### Protein-RNA Correlation (pipeline/modules/cite/cite_corr.py)
- Within-cell-type Spearman r + Fisher z-transform aggregation (June 2026)
- celltype_key="auto" resolves: adt_celltype_manual → adt_celltype →
  cell_type_vote → cell_type
- results_per_celltype DataFrame in return dict
- uns['omicsage_cite_corr_per_ct'] written

### Epitope Characterisation (pipeline/modules/cite/cite_epitope.py)
- Protein-level cell type marker summary
- Epitope coverage table
- Parallel plots show top positive-r pairs (June 2026 fix)

### CITE Pipeline Runner (run_cite_pipeline.py)
- cite_01 through cite_10 step keys
- Per-step HTML reports, cite_00_combined_report.html

### CITE Reports (reports/templates/cite/)
- cite_normalize_report.py, cite_doublets_report.py,
  cite_reduce_report.py, cite_harmony_report.py,
  cite_annotate_report.py (large-cluster flag, confidence table),
  cite_integration_report.py (scib interpretation column,
  preferred embedding verdict, ? batches fix),
  cite_deg_report.py (platelet flag, all-negative flag, ? cells fix),
  cite_gsea_report.py (ribosomal artifact flag, KEGG disease flag),
  cite_corr_report.py (per-cell-type heatmap),
  cite_epitope_report.py
- cite_combined_report.py → cite_00_combined_report.html

### Manual Review (CITE_MASTER_PROMPT.md)
- 10-task structured review prompt for cite_00_combined_report.html
- Attach + send, no pipeline required

### Benchmark Results (GSE194122 BMMC CITE-seq, 21,778 cells)
- 134-protein ADT panel, 8 cell types annotated
- MOFA+ iLISI=0.054, cLISI=0.9997, ASW=0.618
- totalVI iLISI=0.229, cLISI=1.000, ASW=0.712
- Recommended embedding: totalVI (better batch mixing + bio preservation)
- Within-cell-type protein-RNA correlations: median r confirmed negative
  for costimulatory T cell markers (CD5, CD28, CD4) — biologically expected
  (activation-induced internalization in BMMC)

### Config (config/runs/GSE194122_cite.yaml)
- Full CITE pipeline config with annotation_map

---

## ✅ Dependency Files (June 2026 — synced)

### environment.yml
- Added: pydeseq2>=0.4, mofapy2>=0.7, singler>=0.5.0,
  celldex>=0.3.0, beautifulsoup4>=4.12
- scikit-misc pinned to >=0.3

### requirements-ci.txt
- Added: python-pptx==1.0.2, beautifulsoup4>=4.12, matplotlib>=3.8
- Aligned: gseapy==1.1.3, muon>=0.1.6, scvi-tools>=1.1

---

## 🔄 Phase 5 — Multiome (RNA + ATAC jointly) — IN PROGRESS

### Session 1 (next)
- Goal 1: Run RNA pipeline on GSE194122 multiome data
- Goal 2: Build atac_qc.py
- Config: config/runs/GSE194122_multiome.yaml (needs update)
- See NEXT_SESSION.md for full plan

### Planned modules
- atac_qc.py — fragment-based QC metrics (TSS enrichment, nucleosome signal)
- atac_reduce.py — TF-IDF → LSI → UMAP
- atac_annotate.py — gene activity scores + RNA-guided annotation
- multiome_integration.py — MultiVI (scvi-tools) joint embedding
- multiome_deg.py — joint RNA DEG + differential accessibility
- run_multiome_pipeline.py + combined report

### Key architectural decisions
- Integration method: MultiVI (WNN deferred — pynndescent hang unresolved)
- ATAC embedding key: obsm['X_lsi'] (LSI, not PCA)
- Joint embedding key: obsm['X_multivi']
- ATAC annotation key: obs['atac_celltype']
- Batch key: batch (same as CITE, same donors)
- Fragment files: check availability before designing TSS module

---

## What Is NOT Built Yet
- scATAC standalone pipeline (absorbed into Multiome phase)
- Spatial transcriptomics (Phase 6)
- Streamlit UI (Phase 7)
- CLI (Phase 7)
- WNN joint embedding (deferred — pynndescent hang)
- DSB normalization (deferred — no empty droplets in GSE194122)

---

## Processed Data Files
### CITE-seq (GSE194122)
- data/processed/GSE194122/01_qc.h5mu
- data/processed/GSE194122/cite_01_normalize_adt.h5mu
- data/processed/GSE194122/cite_02_doublets.h5mu
- data/processed/GSE194122/cite_03_reduce_adt.h5mu
- data/processed/GSE194122/cite_04_harmony_adt.h5mu
- data/processed/GSE194122/cite_05_annotate_adt.h5mu
- data/processed/GSE194122/cite_06_integration.h5mu
- data/processed/GSE194122/cite_07_deg.h5mu
- data/processed/GSE194122/cite_08_gsea.h5mu
- data/processed/GSE194122/cite_09_corr.h5mu
- data/processed/GSE194122/cite_10_epitope.h5mu

### scRNA (GSE194122 multiome RNA arm)
- data/processed/GSE194122_multiome/ ← to be created in Phase 5 Session 1

### scRNA (GSE166635 HCC benchmark)
- data/processed/GSE166635/05_annotated.h5ad
