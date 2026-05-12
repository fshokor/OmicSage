# OmicSage — Phase 3 AI Layer: Complete Plan
> Written: May 2026  
> Status: Approved design — ready to implement session by session  
> Hard rule: one session = one item on this plan. No exceptions.

---

## The Core Idea

The pipeline runs perfectly without AI. When an expert bioinformatician is available, they review the output themselves. When they are not, the AI layer steps in to do what that expert would do — at every stage of the analysis.

AI in OmicSage is not a chatbot bolted onto a pipeline. It is a structured scientific collaborator that operates in three modes:

| Mode | What it does | When it runs |
|------|-------------|--------------|
| **Decision support** | Looks at your data and tells you what to do next | Before each analysis step |
| **Interpretation** | Looks at your results and explains what they mean biologically | After each analysis step |
| **Communication** | Assembles everything into a deliverable a biologist can act on | At the end of the full run |

All three modes are fully optional and independently controllable. The AI layer has three levels of optionality — the pipeline and reports run perfectly at every level.

| Level | How to control | Effect |
|-------|---------------|--------|
| **Global off** | `ai_features: false` | Nothing in `ai/` is imported. Pipeline and reports unchanged. |
| **Per-module off** | `modules: { cluster_annotator: false }` | That module is skipped. Others still run. Report omits that section. |
| **Per-step at runtime** | `--ai` flag on `run_pipeline.py` | AI runs only for the requested step. |

A module that is disabled or fails always returns `None` silently. The pipeline and report handle `None` gracefully — that section is simply absent from the output. No errors, no broken reports.

---

## Foundation: The Study Context File

Before any AI feature runs, the scientist provides a `study_context.yaml` at the start of the project. This is the single most important input to the AI layer — the more context provided, the better every downstream AI output will be.

### File location
```
config/study_context.yaml   ← filled once per project, referenced by all AI modules
```

### Schema

```yaml
# study_context.yaml
# Fill in what you know. Leave fields blank or set to null if unknown.
# The AI will work with whatever is provided — more context = better output.

dataset:
  name: "GSE166635"
  species: "human"                    # human | mouse | other
  tissue: "liver"                     # required — used in every AI prompt
  cell_source: "tumour + adjacent"    # how cells were collected

disease:
  context: "hepatocellular carcinoma" # null if healthy tissue / not applicable
  stage: null                         # clinical stage if known
  treatment: null                     # treatment history if known

experiment:
  design: "tumour_vs_normal"          # tumour_vs_normal | time_course | perturbation | atlas | null
  n_donors: 6
  n_conditions: 2
  conditions: ["tumour", "adjacent_normal"]
  batch_key: "sample_id"              # obs column used for batch, null if single batch

biological_question: |
  Characterise the immune and stromal microenvironment of hepatocellular carcinoma.
  Identify immunosuppressive cell populations and potential therapeutic targets.
  # Leave blank and the AI will generate a question based on the data type and tissue.

objectives:
  - "Identify major cell types in the tumour microenvironment"
  - "Find differentially expressed genes between tumour and normal conditions"
  - "Highlight potential drug targets"
  # These become the Objectives section of the final report.

notes: |
  Wang et al. 2025 (npj Precision Oncology) published analysis of this dataset.
  Key finding: immunosuppressive macrophage subpopulation driving immune evasion.
  # Any prior knowledge about this dataset the AI should be aware of.
```

### Rules
- `tissue` is the only required field. Everything else is optional.
- If `biological_question` is blank, the AI generates one based on tissue + disease + experiment design before the first analysis step.
- If `design` is null, the pipeline advisor infers it from the data (n_conditions, batch structure).
- The file is human-readable and human-editable. A biologist can fill it in without bioinformatics training.
- A CLI wizard (`omicsage init`) will eventually guide the user through filling this in interactively — build the YAML schema first, wrap the wizard around it later.

---

## Foundation: The Prompt Skills System

Every AI feature in OmicSage uses a **skill** — a structured, versioned, testable prompt template stored in `ai/skills/`.

### Why skills instead of inline prompts

- Prompts are versioned independently of code — improve a prompt without touching pipeline logic
- Skills are testable in isolation — you can evaluate a skill's output quality without running the full pipeline
- Skills are provider-agnostic — the same skill runs on Claude, GPT-4o, or Ollama
- Skills are documented — every skill declares its inputs, outputs, and what it has been tested on

### Skill file format

```yaml
# ai/skills/cluster_annotator.yaml

name: cluster_annotator
version: "1.0"
description: >
  Annotates single-cell clusters by mapping marker genes to known cell types.
  Returns cell type label, confidence, supporting markers, and recommended
  reference database for automated annotation.

tested_on:
  - PBMC (10x Genomics)
  - HCC tumour microenvironment (GSE166635)

inputs:
  - name: tissue
    type: str
    description: Tissue of origin (e.g. "liver", "PBMC", "lung")
  - name: disease_context
    type: str | null
    description: Disease context if applicable
  - name: cluster_id
    type: str
    description: Cluster identifier
  - name: marker_genes
    type: list[str]
    description: Top N marker genes for this cluster, ranked by log2FC
  - name: n_cells
    type: int
    description: Number of cells in this cluster

output_schema:
  cell_type: str
  confidence: "high | medium | low"
  supporting_markers: list[str]   # subset of input markers that drove the call
  alternative_types: list[str]    # other possibilities if ambiguous
  recommended_db: str             # CellTypist | PanglaoDB | Allen Brain | other
  manual_marker_set: list[str]    # genes a biologist would verify manually
  reasoning: str                  # one paragraph, cites markers explicitly

system_prompt: |
  You are an expert computational biologist specialising in single-cell RNA-seq
  analysis and cell type annotation. You have deep knowledge of cell type marker
  genes across tissues and diseases.

  You always ground your annotations in the specific marker genes provided.
  You never invent markers not present in the input list.
  You always recommend the most appropriate reference database for the tissue type.
  You respond only with valid JSON matching the output schema.

user_prompt_template: |
  Annotate the following single-cell cluster.

  Tissue: {tissue}
  Disease context: {disease_context}
  Cluster ID: {cluster_id}
  Number of cells: {n_cells}
  Top marker genes (ranked by log2FC): {marker_genes}

  Return a JSON object with these exact keys:
  - cell_type: your best cell type prediction
  - confidence: high | medium | low
  - supporting_markers: which of the input markers drove your call
  - alternative_types: other possibilities if this is ambiguous (empty list if confident)
  - recommended_db: the best reference database for automated annotation of this tissue
  - manual_marker_set: 5-8 genes a biologist should check manually to verify this call
  - reasoning: one paragraph citing specific markers and your reasoning
```

### Skill loader (Python)

```python
# ai/_skill_loader.py
# Loads a skill YAML, fills the template, returns (system_prompt, user_prompt)
def load_skill(skill_name: str, **kwargs) -> tuple[str, str]:
    ...
```

All feature modules call `load_skill()` — they never contain raw prompt strings.

---

## Shared Infrastructure (Session 0 — build before any feature)

### Files to build

```
ai/
├── _llm_client.py       # BioChatter wrapper, routes on provider: claude|ollama|openai
├── _audit_log.py        # write_audit_record() → logs/llm/<module>.jsonl
├── _skill_loader.py     # loads + fills skill YAML templates
├── _config_gate.py      # raises AiDisabledError if ai_features: false
├── _base.py             # AiResult base dataclass (timestamp, model, skill, reasoning)
└── skills/              # all skill YAML files live here
    ├── pipeline_advisor.yaml
    ├── clustering_advisor.yaml
    ├── cluster_annotator.yaml
    ├── deg_validator.yaml
    ├── coherence_reviewer.yaml
    ├── downstream_suggester.yaml
    ├── narrative_generator.yaml
    └── report_writer.yaml
```

### LLM client routing

```yaml
# config.yaml — full ai section schema
ai:
  features: true              # false = disable everything, no imports, no API key needed
  provider: claude            # claude | ollama | openai
  model: claude-sonnet-4-20250514
  # model: llama3             # for ollama (free, local, no key)
  # model: gpt-4o             # for openai

  # Per-module control — omit a key to default to true
  # Set any module to false to skip it while keeping others active
  modules:
    pipeline_advisor:     true
    clustering_advisor:   true
    cluster_annotator:    true
    deg_validator:        true
    coherence_reviewer:   true
    downstream_suggester: true
    narrative_generator:  true
    report_writer:        true
```

Provider routing happens once in `_llm_client.py`. Every feature module calls `call_llm(skill, inputs)` and never touches BioChatter directly.

### The config gate — three-level check

`_config_gate.py` exposes a single function that every feature module calls as its first line:

```python
from ai._config_gate import check_ai_enabled
check_ai_enabled(config, module="cluster_annotator")
# Raises AiDisabledError if:
#   - config.ai.features is false, OR
#   - config.ai.modules.cluster_annotator is false
# Returns silently if the module is enabled.
```

Callers catch `AiDisabledError` and return `None`. They never catch any other exception — real errors propagate normally so they can be fixed.

```python
# Pattern every feature module follows — no exceptions
def run(adata, config):
    try:
        check_ai_enabled(config, module="cluster_annotator")
    except AiDisabledError:
        return None
    # ... rest of the module
```

### Per-step runtime control

The `--ai` flag on `run_pipeline.py` activates AI only for the requested step:

```bash
# Run full pipeline, no AI
python run_pipeline.py --config config/runs/GSE166635.yaml --step all

# Run only clustering step with AI active
python run_pipeline.py --config config/runs/GSE166635.yaml --step cluster --ai

# Run annotation step with AI, using a specific module override
python run_pipeline.py --config config/runs/GSE166635.yaml --step annotate --ai --ai-modules cluster_annotator
```

Without `--ai`, the pipeline ignores the `ai` section entirely regardless of what is in the config file.

### Audit log format (one line per LLM call)

```jsonl
{
  "timestamp": "2026-05-12T14:32:01Z",
  "module": "cluster_annotator",
  "skill_version": "1.0",
  "model": "claude-sonnet-4-20250514",
  "provider": "claude",
  "input_summary": {"tissue": "liver", "cluster_id": "3", "n_markers": 20},
  "prompt_tokens": 412,
  "completion_tokens": 287,
  "raw_response": "...",
  "parsed_output": {...},
  "parse_success": true
}
```

### BioChatter API reference (verified on v0.14.2)

**Pinned version:** `biochatter==0.14.2` — do not upgrade without re-running all AI tests.
BioChatter has broken its API across minor versions. Two confirmed breaking changes:
- `set_system_message()` was renamed to `append_system_message()`
- `base_url` became a required positional argument for `OllamaConversation`

**Correct usage pattern for every provider:**

```python
# Ollama (local, free)
from biochatter.llm_connect import OllamaConversation
conv = OllamaConversation(
    base_url='http://localhost:11434',
    prompts={},
    model_name='llama3',
    correct=False,
)
conv.append_system_message('You are an expert bioinformatician.')
response, token_usage, correction = conv.query('Your prompt here')

# Claude (Anthropic)
from biochatter.llm_connect import AnthropicConversation
conv = AnthropicConversation(
    model_name='claude-sonnet-4-20250514',
    prompts={},
    correct=False,
)
conv.set_api_key(api_key=os.environ['ANTHROPIC_API_KEY'])
conv.append_system_message('You are an expert bioinformatician.')
response, token_usage, correction = conv.query('Your prompt here')

# OpenAI
from biochatter.llm_connect import GptConversation
conv = GptConversation(
    model_name='gpt-4o',
    prompts={},
    correct=False,
)
conv.set_api_key(api_key=os.environ['OPENAI_API_KEY'])
conv.append_system_message('You are an expert bioinformatician.')
response, token_usage, correction = conv.query('Your prompt here')
```

**Methods confirmed present on v0.14.2** (from `dir(OllamaConversation)`):
`append_system_message`, `append_user_message`, `append_ai_message`,
`query`, `set_api_key`, `set_prompts`, `reset`, `setup`

**`_llm_client.py` must use these exact method names.** Any upgrade to biochatter
requires re-verifying this list before touching any AI module.

### Tests required before any feature session starts
- Config gate raises `AiDisabledError` when `ai_features: false`
- Config gate raises `AiDisabledError` when module key is explicitly `false`
- Config gate passes when module key is missing (defaults to enabled)
- Config gate passes when `ai_features: true` and module is enabled
- `--ai` flag absent → config gate raises regardless of config file contents
- LLM client routes to correct BioChatter class per provider
- Unknown provider raises `ValueError` listing valid options
- Skill loader fills template correctly, raises on missing required input
- Audit log writes and appends correctly
- Mock LLM client injectable for all feature tests

---

## Feature Sessions

### Session 1 — A1: Pipeline Advisor
**File:** `ai/pipeline_advisor.py`  
**Skill:** `ai/skills/pipeline_advisor.yaml`  
**When it runs:** After data loading, before any analysis step

**What it does:**  
Reads the study context + loaded data properties, then recommends which pipeline steps to run and in what order. This is the first thing a user without a bioinformatics background sees.

**Inputs to the skill:**
- From `study_context.yaml`: tissue, disease, experiment design, biological question
- From loaded adata/mdata: n_cells, n_genes, modalities present, n_batches, n_donors, n_conditions

**Output:**
```python
@dataclass
class PipelineAdvice:
    recommended_steps: list[StepRecommendation]
    # Each StepRecommendation: step_name, priority (required|recommended|optional), rationale
    inferred_biological_question: str | None  # if biological_question was blank
    warnings: list[str]                       # e.g. "only 2 samples — pseudobulk unreliable"
    reasoning: str
```

**Key logic (not AI — rule-based first, AI for rationale):**
- n_batches > 1 → recommend batch correction (Harmony or scVI based on n_cells)
- n_donors > 2 AND n_conditions > 1 → recommend pseudobulk DEG over Wilcoxon
- modalities includes ADT → recommend WNN integration (flag as future phase)
- n_cells < 500 → warn against certain steps

**Tests:**
- Returns None when ai_features=False
- Batch detection fires correctly on mock adata with batch_key
- Pseudobulk recommendation fires when n_donors > 2
- inferred_biological_question generated when blank in study_context
- Audit log written

---

### Session 2 — A2: Clustering Advisor
**File:** `ai/clustering_advisor.py`  
**Skill:** `ai/skills/clustering_advisor.yaml`  
**When it runs:** Before clustering step, after dim. reduction

**What it does:**  
Suggests Leiden resolution and approximate expected cluster count, informed by:
1. Literature: what resolutions have published papers used on comparable tissue/disease?
2. Data: silhouette scores and cluster stability across resolution sweep (already computed by pipeline)

**Inputs:**
- tissue, disease_context from study_context
- resolution_sweep_results: dict of resolution → silhouette_score (from pipeline)
- n_cells, n_highly_variable_genes

**Output:**
```python
@dataclass
class ClusteringAdvice:
    suggested_resolution: float
    resolution_range: tuple[float, float]   # reasonable range to explore
    expected_n_clusters: int
    literature_context: list[PubMedRef]     # comparable studies and their resolutions
    reasoning: str
```

**PubMed RAG query pattern:**
`"Leiden clustering resolution {tissue} {disease} single-cell RNA-seq"`

**Tests:**
- Returns None when ai_features=False
- Resolution suggestion within range of input sweep
- PubMed query constructed correctly from study_context fields
- Empty literature result handled gracefully

---

### Session 3 — B1: Cluster Annotator
**File:** `ai/cluster_annotator.py`  
**Skill:** `ai/skills/cluster_annotator.yaml`  
**When it runs:** After clustering + marker gene computation

**What it does:**  
Maps top N marker genes per cluster to known cell types. Returns three outputs per cluster — not just a label:
1. Automated label (for the report and `adata.obs['ai_cell_type']`)
2. Recommended reference database (for the pipeline config — CellTypist, PanglaoDB, Allen Brain, etc.)
3. Manual marker set (for the biologist to verify at the microscope or in literature)

**Inputs:**
- `adata.uns['rank_genes_groups']` — top 20 markers per cluster
- tissue, disease_context from study_context
- cluster cell counts

**Output written to adata:**
- `adata.obs['ai_cell_type']` — predicted label
- `adata.obs['ai_confidence']` — high | medium | low
- `adata.obs['ai_alt_types']` — comma-separated alternatives if ambiguous

**Report section:** UMAP coloured by `ai_cell_type` alongside SingleR `cell_type_vote`

**Tests:**
- Returns None when ai_features=False
- Correct parsing with 3-cluster mock response
- obs columns written correctly
- Missing rank_genes_groups raises informative error
- Partial parse failure (1 of 3 clusters) does not crash — logs warning, skips that cluster

---

### Session 4 — B2: DEG Validator + Literature Linker
**File:** `ai/deg_validator.py`  
**Skill:** `ai/skills/deg_validator.yaml`  
**When it runs:** After DEG computation

**What it does:**  
Two-part analysis per comparison:
1. **Validation:** Are the top DEGs consistent with known biology of this cell type and context? Flags unexpected genes.
2. **Literature linking:** For each top gene, retrieves PMID + one-sentence context via PubMed RAG.

**Inputs:**
- Top 10 DEGs per cluster/comparison (gene symbol, log2FC, adjusted p-value)
- Cell type annotation per cluster (from B1 or SingleR)
- tissue, disease_context from study_context

**Output:**
```python
@dataclass
class DegValidation:
    expected_genes: list[str]        # consistent with known biology + PMID evidence
    unexpected_genes: list[str]      # novel or potentially artifactual — worth investigating
    literature_links: list[GeneLitRef]  # gene → pmid + title + one-sentence context
    validation_summary: str          # overall assessment paragraph
    discovery_highlights: list[str]  # "Gene X upregulated in cluster Y — not previously reported in this context"
```

**Copyright constraint:** Only PMIDs + titles stored and reported. Never abstract text.

**Tests:**
- Expected/unexpected classification on mock DEG list
- Empty gene list returns empty result gracefully
- PMID deduplication across genes
- Literature table renders in report without abstract text

---

### Session 5 — B3: Coherence Reviewer
**File:** `ai/coherence_reviewer.py`  
**Skill:** `ai/skills/coherence_reviewer.yaml`  
**When it runs:** After all analysis steps complete, before report generation

**What it does:**  
Reads a compressed summary of the full analysis and acts as a senior bioinformatician reviewer. Identifies inconsistencies, surprises, and things worth looking deeper into.

**The key engineering challenge — `analysis_summary.json`:**  
Must compress the full run into ~2000 tokens while preserving all information needed for coherent review. Structure:

```json
{
  "study_context": { "tissue": "...", "disease": "...", "n_cells": 5000, "n_batches": 3 },
  "qc_decisions": { "min_genes": 300, "max_mt_pct": 20, "cells_removed_pct": 12.3 },
  "clustering": { "resolution": 0.5, "n_clusters": 14, "silhouette_score": 0.41 },
  "cell_types": [
    { "cluster": 0, "cell_type": "CD8+ T cell", "confidence": "high", "n_cells": 412 },
    ...
  ],
  "top_degs": [
    { "comparison": "tumour_vs_normal", "cluster": 0, "gene": "PDCD1", "log2fc": 2.1 },
    ...
  ],
  "pathways": [
    { "cluster": 0, "pathway": "T cell exhaustion", "padj": 0.001 },
    ...
  ]
}
```

**Output:**
```python
@dataclass
class CoherenceReview:
    flags: list[CoherenceFlag]
    # Each flag: category (qc|clustering|annotation|deg|pathway), severity (info|warning|critical),
    #            description, suggestion
    sub_clustering_candidates: list[str]  # cluster IDs worth splitting
    rare_cell_candidates: list[str]       # potential rare populations to investigate
    overall_assessment: str
```

**Example flags:**
- "Cluster 4 contains both naive and memory T cell markers — consider sub-clustering"
- "QC removed 35% of cells — higher than typical for this tissue, verify thresholds"
- "PDCD1 upregulation in CD8+ T cells consistent with exhaustion phenotype — confirmed in 12 published HCC studies"

**Tests:**
- analysis_summary.json builds correctly from mock adata
- Flag generation from mock analysis with known inconsistency
- Sub-clustering candidate identified when mixed markers present
- Graceful handling of missing analysis steps (e.g. no GSEA results)

---

### Session 6 — A3: Downstream Analysis Suggester
**File:** `ai/downstream_suggester.py`  
**Skill:** `ai/skills/downstream_suggester.yaml`  
**When it runs:** After coherence review, as a forward-looking step

**What it does:**  
Given the annotated, reviewed analysis, suggests what to do next based on what was found and what the original biological question was.

**Logic:**
- Progenitor + mature cell types detected → suggest trajectory analysis (Slingshot/PAGA)
- Immune + tumour cells present → suggest cell-cell communication (CellChat/LIANA)
- Clinical metadata available → suggest survival analysis
- Rare cell candidate flagged by B3 → suggest sub-clustering workflow
- Multiple time points → suggest pseudotime ordering

**Output:**
- Written to `reports/<dataset>/NEXT_STEPS.md` — human-readable markdown, not just machine output
- Each suggestion: step name, rationale, expected output, relevant tool

**Tests:**
- Trajectory suggestion fires when progenitor/mature pair detected
- Communication suggestion fires when immune + non-immune cells present
- Output is valid markdown with correct structure

---

### Session 7 — C1: Biological Narrative Generator
**File:** `ai/narrative_generator.py`  
**Skill:** `ai/skills/narrative_generator.yaml`  
**When it runs:** After all interpretation features have run

**What it does:**  
Assembles outputs from A1, B1, B2, B3 into four narrative blocks for the report:

1. **QC rationale** — why these thresholds, what was removed, whether it is within normal range for this tissue
2. **Cell type landscape** — what populations were found, in what proportions, what is notable given the tissue and disease
3. **Differential expression** — what changed between conditions, what it means biologically, how it connects to the literature
4. **Interpretation + perspectives** — what the findings suggest, what the open questions are, what experiments would validate the key claims

**Groundedness requirement:**  
Every factual claim must cite a metric value, a gene name, or a PMID. Score = cited factual sentences / total factual sentences. Target ≥ 0.85. Measured by `tests/test_groundedness.py`.

**Graceful degradation:**  
If any upstream AI module output is missing, that narrative block is skipped. The report still generates with the remaining blocks.

**Output:**  
- `reports/<dataset>/ai_narrative.md` — standalone markdown file
- Injected into combined HTML report as a new "AI interpretation" tab

---

### Session 8 — C2: Full Report + PowerPoint Generator
**File:** `ai/report_writer.py`  
**Skill:** `ai/skills/report_writer.yaml`  
**When it runs:** Final step of the pipeline

**What it does:**  
Generates two deliverables from the full analysis + narrative:

**HTML/PDF report sections:**
1. Objectives — from study_context biological_question + objectives fields
2. Data — dataset description, n_cells, modalities, QC summary table
3. Methods — auto-generated from config + software versions (Phase 2 engine)
4. Results — all figures + C1 narrative blocks
5. Conclusion — AI-generated summary of the 3–5 key findings
6. Perspectives — AI-generated: open questions, suggested follow-up experiments, potential clinical relevance

**PowerPoint (python-pptx):**
- Title slide
- Data overview slide (n_cells, modalities, QC metrics)
- UMAP slide (coloured by cell type)
- Cell type proportions slide
- Top DEGs slide (volcano or heatmap)
- Pathway enrichment slide
- Key findings slide
- Conclusions + perspectives slide
- Speaker notes on every slide = corresponding C1 narrative paragraph

**This is the freelance deliverable.** A biologist client receives a PowerPoint they can present to their PI the next day.

---

## Phase 3 Milestone Definition

Phase 3 is complete when ALL of the following are true:

| Check | How verified |
|-------|-------------|
| All prior tests still passing (231+) | `python -m pytest tests/ -v` |
| All 8 AI feature modules have passing tests | `python -m pytest tests/test_ai_*.py -v` |
| Full pipeline runs on GSE166635 with `ai_features: true` | `run_pipeline.py --config ... --step all` |
| Full pipeline runs on GSE166635 with `ai_features: false` | Same command, no AI output, no errors |
| Groundedness score ≥ 0.85 on GSE166635 narrative | `tests/test_groundedness.py` |
| All LLM calls written to `logs/llm/*.jsonl` | Check log files exist and are non-empty |
| Ollama provider works locally with llama3 | Manual test with `provider: ollama` |
| PowerPoint generated with all 8 slides | Open and verify manually |
| `NEXT_STEPS.md` written to run output directory | Check file exists and is non-empty |
| `analysis_summary.json` written per run | Check file exists and is valid JSON |

---

## Build Order Summary

```
Session 0 — Infrastructure: _llm_client, _audit_log, _skill_loader, _config_gate, _base
Session 1 — A1: Pipeline advisor
Session 2 — A2: Clustering advisor (first use of PubMed RAG)
Session 3 — B1: Cluster annotator
Session 4 — B2: DEG validator + literature linker
Session 5 — B3: Coherence reviewer (build analysis_summary.json here)
Session 6 — A3: Downstream analysis suggester
Session 7 — C1: Narrative generator
Session 8 — C2: Full report + PowerPoint
Session 9 — Milestone validation: groundedness test + end-to-end GSE166635
```

Total: 10 sessions. At 1 session per day (your schedule), Phase 3 completes in 2 weeks.

---

## Files Created by Phase 3

```
ai/
├── _llm_client.py
├── _audit_log.py
├── _skill_loader.py
├── _config_gate.py
├── _base.py
├── pipeline_advisor.py
├── clustering_advisor.py
├── cluster_annotator.py
├── deg_validator.py
├── coherence_reviewer.py
├── downstream_suggester.py
├── narrative_generator.py
├── report_writer.py
└── skills/
    ├── pipeline_advisor.yaml
    ├── clustering_advisor.yaml
    ├── cluster_annotator.yaml
    ├── deg_validator.yaml
    ├── coherence_reviewer.yaml
    ├── downstream_suggester.yaml
    ├── narrative_generator.yaml
    └── report_writer.yaml

config/
└── study_context.yaml          ← filled once per project

reports/<dataset>/
├── ai_narrative.md             ← written by C1
├── NEXT_STEPS.md               ← written by A3
└── analysis_summary.json       ← written by B3, consumed by C1 + C2

logs/llm/
├── pipeline_advisor.jsonl
├── clustering_advisor.jsonl
├── cluster_annotator.jsonl
├── deg_validator.jsonl
├── coherence_reviewer.jsonl
├── downstream_suggester.jsonl
├── narrative_generator.jsonl
└── report_writer.jsonl

tests/
├── test_ai_infrastructure.py
├── test_pipeline_advisor.py
├── test_clustering_advisor.py
├── test_cluster_annotator.py
├── test_deg_validator.py
├── test_coherence_reviewer.py
├── test_downstream_suggester.py
├── test_narrative_generator.py
├── test_report_writer.py
└── test_groundedness.py
```

---

## Design Principles (carry into every session)

1. **Skills are the prompts.** No raw strings in Python files. Every prompt lives in `ai/skills/` as a versioned YAML.
2. **AI is optional at three levels.** Global (`ai_features: false`), per-module (`modules: { x: false }`), and per-step (`--ai` flag at runtime). A disabled module always returns `None` silently. The pipeline never requires AI to run.
3. **Every LLM call is audited.** No exceptions. Timestamp, model, inputs, raw response, parse result.
4. **Graceful degradation at every level.** A failed AI module writes a warning to the log and returns `None`. It never crashes the pipeline.
5. **Context scales quality.** The more the scientist fills in `study_context.yaml`, the better every AI output. The AI never blocks on missing context.
6. **The analysis_summary.json is load-bearing.** B3 and C1 both depend on it. Design it carefully — it must compress the full run into ~2000 tokens without losing the information needed for coherent review.
7. **Copyright is enforced at the report layer.** PubMed results provide PMIDs and titles only. Never abstract text in any output file.
