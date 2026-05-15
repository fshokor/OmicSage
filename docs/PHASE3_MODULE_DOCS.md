# OmicSage — Phase 3 AI Layer: Module Documentation
> Created: 2026-05-15
> Updated: after each Phase 3 session
> Scope: ai/ directory only — infrastructure + all feature modules
> Companion to: .dev_memory/PHASE3_PLAN.md (design decisions + session specs)

---

## How the AI Layer Works

The AI layer sits on top of the complete Phase 1-2 pipeline and is fully
optional. Setting `ai_features: false` in config.yaml disables it entirely —
the pipeline and reports run unchanged with zero API calls.

Every feature module follows the same four-step pattern:

```
1. check_ai_enabled()     ← gate check — returns None if disabled
2. rule-based pre-checks  ← fast logic that runs before any LLM call
3. call_llm()             ← loads skill YAML, calls BioChatter, writes audit log
4. parse + return         ← JSON parse → typed dataclass inheriting AiResult
```

The skill YAML is the prompt. Python files contain zero raw strings.

---

## Infrastructure Modules (Session 0 — do not modify)

### ai/_base.py
**Status:** ✅ Complete — 3 tests passing

Base dataclass for all AI feature outputs.

```python
@dataclass
class AiResult:
    timestamp: str    # ISO-8601 UTC, auto-populated at instantiation
    model: str        # e.g. "claude-sonnet-4-20250514"
    provider: str     # "claude" | "ollama" | "openai"
    skill_name: str   # e.g. "cluster_annotator"
    skill_version: str
    reasoning: str    # free-text reasoning from LLM response
```

All feature dataclasses inherit from AiResult. Subclass example:

```python
@dataclass
class PipelineAdvice(AiResult):
    recommended_steps: list[StepRecommendation] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
```

---

### ai/_config_gate.py
**Status:** ✅ Complete — 7 tests passing

Three-level optionality check. Called as the first line of every feature module.

```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config, *, runtime_ai=True):
    try:
        check_ai_enabled(config, module="my_module", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None   # silent — no log, no warning
```

**Three levels checked in order:**

| Level | What triggers it | Config location |
|-------|-----------------|-----------------|
| 1 — global | `ai.features: false` | config.yaml top-level |
| 2 — per-module | `ai.modules.my_module: false` | config.yaml ai.modules dict |
| 3 — runtime | `--ai` flag absent on CLI | `runtime_ai=False` passed by runner |

Missing module key → defaults to **enabled** (opt-out model, not opt-in).

**Public API:**
```python
check_ai_enabled(config: dict, module: str, runtime_ai: bool = True) -> None
# Returns None silently or raises AiDisabledError
```

---

### ai/_audit_log.py
**Status:** ✅ Complete — 4 tests passing

Append-only JSONL audit log. One file per module under `logs/llm/`.

```python
write_audit_record(
    log_dir="logs/llm",
    module="cluster_annotator",
    skill_version="1.0",
    model="claude-sonnet-4-20250514",
    provider="claude",
    input_summary={"tissue": "liver", "cluster_id": "3"},
    token_usage={"prompt_tokens": 412, "completion_tokens": 287},
    raw_response=raw_str,
    parsed_output=result_dict,   # None if parse failed
    parse_success=True,
)
```

**Record schema (one JSON object per line):**
```
timestamp          ISO-8601 UTC string
module             str
skill_version      str
model              str
provider           str
input_summary      dict   (scalar values only — no AnnData objects)
prompt_tokens      int | null
completion_tokens  int | null
raw_response       str    (full LLM output before any parsing)
parsed_output      dict | null
parse_success      bool
```

**Never raises.** A write failure prints to stderr and continues.
Log dir is created automatically if absent.

---

### ai/_skill_loader.py
**Status:** ✅ Complete (built prior session) — all checks passing

Loads a skill YAML and fills the user prompt template.

```python
from ai._skill_loader import load_skill

system_prompt, user_prompt = load_skill(
    "cluster_annotator",
    tissue="liver",
    cluster_id="3",
    marker_genes=["PDCD1", "CD8A", "GZMB"],
    n_cells=412,
    disease_context="hepatocellular carcinoma",
)
```

- Looks up `ai/skills/<skill_name>.yaml`
- Returns `(system_prompt: str, user_prompt: str)` — ready to pass to BioChatter
- Raises `KeyError` if a required template variable is missing from inputs

---

### ai/_llm_client.py
**Status:** ✅ Complete — 7 tests passing

BioChatter wrapper. The only file that imports BioChatter directly.

**Public API:**
```python
raw_response: str = call_llm(
    skill_name="cluster_annotator",
    inputs={"tissue": "liver", ...},
    config=config,
    log_dir="logs/llm",        # default
    module="cluster_annotator", # defaults to skill_name
    runtime_ai=True,            # default
)
```

**Provider routing:**
```
config.ai.provider = "claude"  → AnthropicConversation (needs ANTHROPIC_API_KEY)
config.ai.provider = "ollama"  → OllamaConversation    (needs local Ollama server)
config.ai.provider = "openai"  → GptConversation       (needs OPENAI_API_KEY)
```

Unknown provider → `ValueError` listing valid options.

**BioChatter version:** `==0.14.2` — pinned. Do not upgrade.
Verified method names on v0.14.2:
- `append_system_message(text)` — NOT `set_system_message` (renamed in v0.14)
- `query(text)` → `(response, token_usage, correction)`
- `set_api_key(api_key=...)` — Claude and OpenAI only

`_build_conversation()` is exposed for monkeypatching in tests.
Feature modules never call BioChatter directly — always via `call_llm()`.

---

### ai/skills/cluster_annotator.yaml
**Status:** ✅ Complete (reference pattern — also used by Session 3)

Establishes the YAML schema that all skill files follow.

**Required top-level keys:**
```yaml
name: str
version: str          # "1.0" — increment when prompt changes
description: str
tested_on: list[str]  # dataset names this has been validated on
inputs:               # list of {name, type, description}
output_schema:        # dict — keys the LLM must return
system_prompt: str    # static — no template variables
user_prompt_template: str  # uses {variable} placeholders
```

---

### config/study_context_template.yaml
**Status:** ✅ Complete

Filled once per project. Referenced by all AI modules.
Copy to `config/runs/<dataset>/study_context.yaml` and fill in.

Required field: `tissue`
All other fields optional — AI works with whatever is provided.

Key fields:
```yaml
dataset:    name, species, tissue, cell_source
disease:    context, stage, treatment
experiment: design, n_donors, n_conditions, conditions, batch_key
biological_question: str   # auto-generated if blank
objectives: list[str]
notes: str
```

---

## Feature Modules

### ai/pipeline_advisor.py — Session 1
**Status:** ✅ Complete

**When it runs:** After data loading, before any analysis step.

**What it does:** Reads study context + loaded data properties, recommends
which pipeline steps to run and in what order. First thing a user without
bioinformatics training sees.

**Output dataclass:**
```python
@dataclass
class StepRecommendation:
    step_name: str
    priority: str   # "required" | "recommended" | "optional"
    rationale: str

@dataclass
class PipelineAdvice(AiResult):
    recommended_steps: list[StepRecommendation]
    inferred_biological_question: str | None   # generated if blank in study_context
    warnings: list[str]
```

**Rule-based pre-checks (before any LLM call):**
- `n_batches > 1` → warn if batch_key not set
- `n_donors > 2` AND `n_conditions > 1` → recommend pseudobulk over Wilcoxon
- modalities includes "ADT" → recommend WNN (flag as Phase 6)
- `n_cells < 500` → warn about unreliable clustering

**Skill file:** `ai/skills/pipeline_advisor.yaml`
**Test file:** `tests/test_pipeline_advisor.py` (12 tests, no real API key)

---

### ai/clustering_advisor.py — Session 2
**Status:** 🔲 Not yet built

**When it runs:** Before clustering, after dimensionality reduction.

**What it does:** Suggests Leiden resolution and expected cluster count,
informed by silhouette scores from the resolution sweep and PubMed RAG
(comparable published studies).

**First use of PubMed RAG** — query pattern:
`"Leiden clustering resolution {tissue} {disease} single-cell RNA-seq"`

**Output dataclass:**
```python
@dataclass
class ClusteringAdvice(AiResult):
    suggested_resolution: float
    resolution_range: tuple[float, float]
    expected_n_clusters: int
    literature_context: list[PubMedRef]   # {pmid, title, context_sentence}
    reasoning: str
```

**Skill file:** `ai/skills/clustering_advisor.yaml`
**Test file:** `tests/test_clustering_advisor.py`

---

### ai/cluster_annotator.py — Session 3
**Status:** 🔲 Not yet built

**When it runs:** After clustering + marker gene computation.

**What it does:** Maps top N marker genes per cluster to known cell types.
Returns label, confidence, recommended reference DB, and manual marker set
for the biologist to verify.

**Writes to adata:**
- `obs['ai_cell_type']` — predicted label
- `obs['ai_confidence']` — "high" | "medium" | "low"
- `obs['ai_alt_types']` — comma-separated alternatives if ambiguous

**Output dataclass (per cluster):**
```python
@dataclass
class ClusterAnnotation(AiResult):
    cluster_id: str
    cell_type: str
    confidence: str
    supporting_markers: list[str]
    alternative_types: list[str]
    recommended_db: str          # CellTypist | PanglaoDB | Allen Brain | other
    manual_marker_set: list[str] # 5-8 genes for biologist to verify
```

**Skill file:** `ai/skills/cluster_annotator.yaml` (already drafted — refine here)
**Test file:** `tests/test_cluster_annotator.py`

---

### ai/deg_validator.py — Session 4
**Status:** 🔲 Not yet built

**When it runs:** After DEG computation.

**What it does:** Two-part analysis per comparison:
1. Validation — are top DEGs consistent with known biology?
2. Literature linking — PMID + one-sentence context per top gene via PubMed RAG.

**Copyright constraint:** PMIDs + titles only. Never abstract text in any output.

**Output dataclass:**
```python
@dataclass
class DegValidation(AiResult):
    expected_genes: list[str]           # consistent with known biology
    unexpected_genes: list[str]         # novel or potentially artifactual
    literature_links: list[GeneLitRef]  # {gene, pmid, title, context_sentence}
    validation_summary: str
    discovery_highlights: list[str]
```

**Skill file:** `ai/skills/deg_validator.yaml`
**Test file:** `tests/test_deg_validator.py`

---

### ai/coherence_reviewer.py — Session 5
**Status:** 🔲 Not yet built

**When it runs:** After all analysis steps, before report generation.

**What it does:** Acts as a senior bioinformatician reviewer. Reads a
compressed `analysis_summary.json` (~2000 tokens) and flags inconsistencies,
surprises, and sub-clustering candidates.

**Load-bearing artifact — analysis_summary.json:**
Built in this session. Consumed by C1 (narrative_generator) and C2 (report_writer).
Schema designed to compress the full run without losing review-relevant information.

```json
{
  "study_context": {"tissue": "...", "n_cells": 5000, ...},
  "qc_decisions": {"min_genes": 300, "cells_removed_pct": 12.3, ...},
  "clustering": {"resolution": 0.5, "n_clusters": 14, ...},
  "cell_types": [{"cluster": 0, "cell_type": "CD8+ T cell", ...}],
  "top_degs": [{"comparison": "...", "gene": "PDCD1", "log2fc": 2.1}],
  "pathways": [{"cluster": 0, "pathway": "T cell exhaustion", ...}]
}
```

**Output dataclass:**
```python
@dataclass
class CoherenceFlag:
    category: str    # "qc" | "clustering" | "annotation" | "deg" | "pathway"
    severity: str    # "info" | "warning" | "critical"
    description: str
    suggestion: str

@dataclass
class CoherenceReview(AiResult):
    flags: list[CoherenceFlag]
    sub_clustering_candidates: list[str]   # cluster IDs worth splitting
    rare_cell_candidates: list[str]
    overall_assessment: str
```

**Skill file:** `ai/skills/coherence_reviewer.yaml`
**Test file:** `tests/test_coherence_reviewer.py`

---

### ai/downstream_suggester.py — Session 6
**Status:** 🔲 Not yet built

**When it runs:** After coherence review, as forward-looking step.

**What it does:** Given the annotated analysis and original biological question,
suggests what to do next. Written to `reports/<dataset>/NEXT_STEPS.md`.

**Trigger logic (rule-based before LLM):**
- Progenitor + mature cell types → suggest trajectory (Slingshot/PAGA)
- Immune + tumour cells → suggest cell-cell communication (CellChat/LIANA)
- Clinical metadata available → suggest survival analysis
- Rare cell candidate from coherence review → suggest sub-clustering
- Multiple time points → suggest pseudotime ordering

**Output:** `reports/<dataset>/NEXT_STEPS.md` — markdown, human-readable.
Each suggestion: step name, rationale, expected output, tool.

**Skill file:** `ai/skills/downstream_suggester.yaml`
**Test file:** `tests/test_downstream_suggester.py`

---

### ai/narrative_generator.py — Session 7
**Status:** 🔲 Not yet built

**When it runs:** After all interpretation modules have run.

**What it does:** Assembles outputs from A1, B1, B2, B3 into four narrative
blocks injected into the combined HTML report as an "AI interpretation" tab.

**Four blocks:**
1. QC rationale — threshold choices, what was removed, normal range for tissue
2. Cell type landscape — populations found, proportions, what is notable
3. Differential expression — what changed, biological meaning, literature links
4. Interpretation + perspectives — key findings, open questions, validation experiments

**Groundedness requirement:**
Every factual claim must cite a metric, gene name, or PMID.
Score = cited factual sentences / total factual sentences. Target ≥ 0.85.
Measured by `tests/test_groundedness.py` (built in Session 9).

**Graceful degradation:** Missing upstream module output → that block is skipped.
Report generates with remaining blocks.

**Output files:**
- `reports/<dataset>/ai_narrative.md`
- Injected into `00_combined_report.html` as new "AI interpretation" tab

**Skill file:** `ai/skills/narrative_generator.yaml`
**Test file:** `tests/test_narrative_generator.py`

---

### ai/report_writer.py — Session 8
**Status:** 🔲 Not yet built

**When it runs:** Final step of the pipeline.

**What it does:** Generates the two client-facing deliverables from the full
analysis + narrative. This is the freelance deliverable — a biologist client
receives a PowerPoint they can present to their PI the next day.

**HTML/PDF report sections:**
1. Objectives — from study_context biological_question + objectives
2. Data — n_cells, modalities, QC summary table
3. Methods — auto-generated from config + software versions
4. Results — all figures + narrative blocks from C1
5. Conclusion — AI summary of 3-5 key findings
6. Perspectives — open questions, follow-up experiments, clinical relevance

**PowerPoint (python-pptx) — 8 slides:**
1. Title
2. Data overview (n_cells, modalities, QC metrics)
3. UMAP coloured by cell type
4. Cell type proportions
5. Top DEGs (volcano or heatmap)
6. Pathway enrichment
7. Key findings
8. Conclusions + perspectives

Speaker notes on every slide = corresponding narrative paragraph from C1.

**Skill file:** `ai/skills/report_writer.yaml`
**Test file:** `tests/test_report_writer.py`

---

## Output Files Written by Phase 3

```
logs/llm/
├── pipeline_advisor.jsonl
├── clustering_advisor.jsonl
├── cluster_annotator.jsonl
├── deg_validator.jsonl
├── coherence_reviewer.jsonl
├── downstream_suggester.jsonl
├── narrative_generator.jsonl
└── report_writer.jsonl

reports/<dataset>/
├── ai_narrative.md          ← Session 7 (C1)
├── NEXT_STEPS.md            ← Session 6 (A3)
└── analysis_summary.json    ← Session 5 (B3), consumed by C1 + C2
```

---

## Config Reference

```yaml
ai:
  features: true              # false = disable everything, no imports, no API key
  provider: claude            # claude | ollama | openai
  model: claude-sonnet-4-20250514
  # model: llama3             # for ollama (free, local, no key)
  # model: gpt-4o             # for openai

  modules:                    # omit a key = defaults to true (opt-out model)
    pipeline_advisor:     true
    clustering_advisor:   true
    cluster_annotator:    true
    deg_validator:        true
    coherence_reviewer:   true
    downstream_suggester: true
    narrative_generator:  true
    report_writer:        true
```

Environment variables:
```
ANTHROPIC_API_KEY   required for provider: claude
OPENAI_API_KEY      required for provider: openai
(Ollama needs no key — connects to http://localhost:11434 by default)
```

---

## Test Count by Module

| Module | Test file | Tests |
|--------|-----------|-------|
| _base + _config_gate + _audit_log + _llm_client | test_ai_infrastructure.py | 20 ✅ |
| pipeline_advisor | test_pipeline_advisor.py | 12 🔲 |
| clustering_advisor | test_clustering_advisor.py | ~10 🔲 |
| cluster_annotator | test_cluster_annotator.py | ~10 🔲 |
| deg_validator | test_deg_validator.py | ~8 🔲 |
| coherence_reviewer | test_coherence_reviewer.py | ~8 🔲 |
| downstream_suggester | test_downstream_suggester.py | ~6 🔲 |
| narrative_generator | test_narrative_generator.py | ~8 🔲 |
| report_writer | test_report_writer.py | ~8 🔲 |
| groundedness | test_groundedness.py | ~4 🔲 |
| **Total Phase 3** | | **~94** |

---

## Design Rules (enforced every session)

1. **Skills are the prompts.** No raw strings in Python. Every prompt lives
   in `ai/skills/` as a versioned YAML.

2. **Three-level optionality.** Global (`ai_features: false`), per-module
   (`modules: { x: false }`), per-step (`--ai` runtime flag). A disabled
   module always returns `None` silently.

3. **Every LLM call is audited.** `write_audit_record()` called after every
   `conv.query()`. No exceptions.

4. **Graceful degradation.** A failed module logs a warning and returns `None`.
   It never crashes the pipeline. Reports omit missing sections gracefully.

5. **Copyright at the report layer.** PubMed: PMIDs + titles only.
   Never abstract text in any output file.

6. **analysis_summary.json is load-bearing.** ~2000 token budget.
   Design it carefully in Session 5 — both C1 and C2 depend on it.

7. **BioChatter pinned at 0.14.2.** Do not upgrade without re-running all
   AI tests and re-verifying method names.
