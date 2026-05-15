# OmicSage — Current Status
> Last updated: 2026-05-15 (Phase 3 Session 0 complete)

## Phase
Phase 3 — AI Layer (Phase 2 complete ✅)

## What Is Built and Tested Right Now

### ✅ Phase 0 — Foundation
- Repo structure scaffolded
- Docker base images (Python + R) — defined, not yet built locally
- CI/CD via GitHub Actions — configured
- YAML config schema — defined
- .dev_memory/ system — initialized

### ✅ Ingestion (pipeline/modules/qc/ingest.py)
- Auto-detects 10x MEX, H5, AnnData formats
- Moves normalized values out of .X, puts raw counts in .X
- Handles GSE194122 CITE-seq and multiome, GSE166635 HCC
- load_dataset_dir() — scans parent folder, loads all MTX subfolders,
  concatenates with obs['sample'] + obs['batch'] per subfolder name
- load_dataset() auto-routes to load_dataset_dir() when given a parent
  folder containing MTX subfolders (no config change needed)

### ✅ QC (pipeline/modules/qc/qc.py)
- Modality-aware: auto-detects GEX / ADT / ATAC from var['feature_types']
- Returns MuData with mdata["rna"], mdata["adt"], mdata["atac"] slots
- MT%, genes/cell, Scrublet doublets
- 42 tests passing in tests/test_qc.py

### ✅ Normalization (pipeline/modules/qc/normalize.py)
- Input: AnnData with raw counts in .X (from mdata["rna"])
- Saves raw counts to layers['counts']
- Normalizes to CP10K, log1p transform
- Saves log1p values to layers['logcounts'] (Seurat convention)
- HVG selection — top 2000 genes, flavor='seurat_v3' (pre-log1p on counts)
- batch_key support for per-batch HVG selection
- Provenance stored in uns['omicsage_normalization']
- 12 tests passing in tests/test_normalize.py

### ✅ Normalization Report (reports/normalization_report.py)
- Self-contained HTML report
- Figures: HVG scatter, library size violin, top 20 HVGs, gene detection rate

### ✅ Dimensionality Reduction (pipeline/modules/qc/reduce.py)
- PCA on HVG subset only (n_comps=50, svd_solver='arpack')
- Auto-select n_pcs via elbow detection (kneed) — falls back to variance threshold
- Neighbor graph (sc.pp.neighbors, n_neighbors=15), UMAP (always), t-SNE (optional)
- Provenance stored in uns['omicsage_reduce']
- 12 tests passing in tests/test_reduce.py

### ✅ Dimensionality Reduction Report (reports/reduce_report.py)
- Figures: scree plot, UMAP × QC metrics (2×2), PCA × QC, optional batch UMAP

### ✅ Leiden Clustering (pipeline/modules/clustering/cluster.py)
- Leiden sweep across configurable resolution_range
- Silhouette score per resolution; best_resolution_override for manual pinning
- obs['leiden_*'] per resolution + obs['leiden'] convenience key
- neighbors_key param — routes to Harmony graph when set
- cluster_key param — stores results in obs['leiden_harmony']
- compute_ari(adata, key_a, key_b) — ARI comparison between any two obs columns
- Provenance stored in uns['omicsage_cluster']
- 16 tests passing in tests/test_cluster.py

### ✅ Clustering Report (reports/cluster_report.py)
- Figures: UMAP grid per resolution, silhouette bar chart, cluster size distribution,
  optional ground-truth cell_type UMAP

### ✅ Cell-Type Annotation (pipeline/modules/annotation/annotate.py)
- Methods: CellTypist (Immune_All_High + Immune_All_Low), marker gene scoring, majority vote
- CellTypist models cached in data/references/celltypist/
- obs columns: celltypist_coarse, celltypist_fine, cell_type_markers,
  cell_type_groundtruth, cell_type_vote, cell_type_confidence
- Provenance stored in uns['omicsage_annotate']
- 18 tests passing, 1 skipped in tests/test_annotate.py

### ✅ Annotation Report (reports/annotate_report.py)
- Figures: UMAP × consensus vote, UMAP × CellTypist fine, confidence distribution
- Per-cluster table with all method labels and confidence scores

### ✅ DEG (pipeline/modules/downstream/deg.py)
- Wilcoxon rank-sum, one-vs-rest, rankby_abs=True, BH correction
- n_genes default 500, exclude_gene_prefixes param
- Provenance stored in uns['omicsage_deg']
- 11 tests passing in tests/test_deg.py

### ✅ DEG Report (reports/deg_report.py)
- Volcano plots + dot plot + summary table + direction column

### ✅ GSEA (pipeline/modules/downstream/gsea.py)
- ORA via gseapy.enrichr — GO BP / KEGG / Reactome
- direction param: "up" | "down" | "both"
- Provenance stored in uns['omicsage_gsea']
- 8 tests passing in tests/test_gsea.py (Enrichr calls mocked — CI-safe)

### ✅ GSEA Report (reports/gsea_report.py)
- Bar charts + bubble plot + direction badges

### ✅ Harmony Batch Correction (pipeline/modules/integration/harmony_correct.py)
- obsm['X_pca_harmony'], obsm['X_umap_precorrection'], obsm['X_umap_harmony']
- neighbors_harmony graph in uns + obsp
- Provenance stored in uns['omicsage_harmony']
- 13 tests passing in tests/test_harmony.py

### ✅ Harmony Report (reports/harmony_report.py)
- Batch composition, mixing metrics, PC shift, UMAP comparison

### ✅ Pseudobulk DEG (pipeline/modules/downstream/pseudobulk_deg.py)
- DESeq2 Wald tests via pydeseq2, one-vs-rest per cell type
- Aggregate layers['counts'] per (cell_type, donor)
- min_cells + min_samples filters with graceful skip + UserWarning
- Output schema identical to deg.py deg_dict
- Provenance stored in uns['omicsage_pseudobulk_deg']
- 14 tests passing in tests/test_pseudobulk_deg.py

### ✅ Pseudobulk DEG Report (reports/pseudobulk_deg_report.py)
- Skipped groups section + pseudobulk-specific stat cards

### ✅ Notebook (notebooks/phase1_qc.ipynb)
- Steps 1–10 complete

### ✅ MILESTONE: Wang et al. 2025 HCC Benchmark
- Full Phase 1 pipeline run on GSE166635 (HCC1 normal + HCC2 tumour)
- Cell types identified: Hepatocytes, T cells, Macrophages, Endothelial,
  Fibroblasts, B cells
- Known HCC markers recovered in DEG results
- Liver metabolism and immune pathway terms in GSEA results

### ✅ Generic Pipeline Runner (run_pipeline.py)
- Config-driven, --from-step / --to-step / --step flags
- Validation at startup, caching, resolution_override support

### ✅ Config System (config/)
- config/schema.yaml — platform master schema
- config/runs/GSE194122.yaml, GSE166635.yaml, GSE194122_multiome.yaml

### ✅ Combined Report (reports/combined_report.py)
- Tabbed HTML assembled from all step reports after every run
- Output: 00_combined_report.html

### ✅ MILESTONE: Phase 2 Complete
- Full pipeline + combined tabbed report from one command on GSE166635

---

## ✅ Phase 3 Session 0 — AI Infrastructure (COMPLETE)

### ai/_base.py
- AiResult base dataclass
- Fields: timestamp (ISO-8601 UTC, auto), model, provider,
  skill_name, skill_version, reasoning
- All feature dataclasses inherit from AiResult
- Tests: 3 passing in test_ai_infrastructure.py

### ai/_config_gate.py
- check_ai_enabled(config, module, runtime_ai=True)
- Raises AiDisabledError at three levels:
    Level 1 — ai.features: false (global)
    Level 2 — ai.modules.<name>: false (per-module)
    Level 3 — runtime_ai=False (--ai flag absent at CLI)
- Missing module key defaults to ENABLED (opt-out model)
- Tests: 7 passing in test_ai_infrastructure.py

### ai/_audit_log.py
- write_audit_record() — appends one JSONL line per LLM call
- Output: logs/llm/<module>.jsonl
- Creates log_dir automatically if absent
- Never raises — logging failure prints to stderr and continues
- Record fields: timestamp, module, skill_version, model, provider,
  input_summary, prompt_tokens, completion_tokens, raw_response,
  parsed_output, parse_success
- Tests: 4 passing in test_ai_infrastructure.py

### ai/_llm_client.py
- call_llm(skill_name, inputs, config, *, log_dir, module, runtime_ai) → str
- _build_conversation(provider, model, config) — injectable for tests
- Routes: claude → AnthropicConversation, ollama → OllamaConversation,
  openai → GptConversation
- Unknown provider raises ValueError listing valid options
- BioChatter version: ==0.14.2 (pinned)
- Verified method names: append_system_message(), query(), set_api_key()
- Tests: 7 passing in test_ai_infrastructure.py

### ai/_skill_loader.py (built prior session)
- load_skill(skill_name, **inputs) → (system_prompt, user_prompt)
- Loads YAML from ai/skills/<name>.yaml
- Fills user_prompt_template with inputs via str.format_map()
- Raises on missing required inputs

### ai/skills/cluster_annotator.yaml (reference pattern)
- Establishes YAML schema all skill files follow
- Fields: name, version, description, tested_on, inputs, output_schema,
  system_prompt, user_prompt_template

### Supporting files
- config/study_context_template.yaml — fill once per project
- config/runs/GSE166635/study_context.yaml — copy + fill
- config/runs/GSE194122/study_context.yaml — copy + fill

### Tests
- tests/test_ai_infrastructure.py — 20 tests, all passing
- test_phase0_structure.py encoding fix: schema fixture now uses encoding="utf-8"

---

## Phase 3 — What Is Being Built Next

Session 1 (next): A1 — Pipeline Advisor
  - ai/skills/pipeline_advisor.yaml
  - ai/pipeline_advisor.py
  - tests/test_pipeline_advisor.py

Remaining sessions:
  Session 2  — A2: Clustering advisor (first PubMed RAG use)
  Session 3  — B1: Cluster annotator
  Session 4  — B2: DEG validator + literature linker
  Session 5  — B3: Coherence reviewer
  Session 6  — A3: Downstream analysis suggester
  Session 7  — C1: Narrative generator
  Session 8  — C2: Full report + PowerPoint
  Session 9  — Milestone validation

## Total Tests Passing
~251 (231 Phase 1-2 + 20 Phase 3 infrastructure)

## What Is NOT Built Yet
- Phase 3: Pipeline advisor (Session 1 — next)
- Phase 3: Clustering advisor
- Phase 3: Cluster annotator
- Phase 3: DEG validator + literature linker
- Phase 3: Coherence reviewer
- Phase 3: Downstream analysis suggester
- Phase 3: Narrative generator
- Phase 3: Full report + PowerPoint
- scVI batch correction → deferred to Phase 6
- ADT QC + CLR normalization
- scATAC module (Phase 4)
- Spatial module (Phase 5)
- Multiome module (Phase 6)
- Streamlit UI (Phase 7)
- CLI (Phase 7)

## Processed Data Files
- data/processed/GSE194122_cite_rna_qc.h5ad
- data/processed/GSE194122_cite_adt_qc.h5ad
- data/processed/GSE194122_multiome_rna_qc.h5ad
- data/processed/GSE194122_multiome_atac_qc.h5ad
- data/processed/GSE166635_HCC1_qc.h5ad
- data/processed/GSE194122_cite_normalized.h5ad
- data/processed/GSE194122_cite_reduced.h5ad
- data/processed/GSE194122_cite_clustered.h5ad
- data/processed/GSE194122_cite_annotated.h5ad
- data/processed/GSE194122_cite_deg.h5ad
- data/processed/GSE194122_cite_gsea.h5ad
- data/processed/GSE194122_cite_harmony.h5ad
- data/processed/GSE194122_cite_harmony_clustered.h5ad
- data/processed/GSE194122_cite_pseudobulk_deg.h5ad
