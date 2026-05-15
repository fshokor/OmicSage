## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 3 — B1 Cluster Annotator built and tested.
                      23 new tests passing (330 total, 1 skipped).
                      Fixed: AiResult base field is skill_name (not skill).
File last worked on: tests/test_cluster_annotator.py

## Today's Goal
Phase 3 — Session 4: B2 — DEG Validator + Literature Linker.

First module that runs PubMed RAG. Nothing else.

Files to create:
  - ai/skills/deg_validator.yaml       ← skill YAML
  - ai/deg_validator.py                ← feature module
  - tests/test_deg_validator.py        ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 330 passed, 1 skipped
```

## Step 2 — Read the DEG validator spec in PHASE3_PLAN.md
Section: "Session 4 — B2: DEG Validator + Literature Linker"
Read it fully before writing any code.

## Step 3 — Build ai/skills/deg_validator.yaml
Follow cluster_annotator.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - tissue, disease_context (from study_context)
  - comparison: str (e.g. "tumour_vs_normal")
  - degs: list of dicts — each with gene, log2fc, padj, cell_type

Output schema (JSON):
  - expected_genes: list of str
  - unexpected_genes: list of str
  - validation_summary: str
  - discovery_highlights: list of str

## Step 4 — Build ai/deg_validator.py

### Public API
```python
from dataclasses import dataclass, field
from ai._base import AiResult

@dataclass
class GeneLitRef:
    gene: str = ""
    pmid: str = ""
    title: str = ""
    context: str = ""   # one-sentence context from PubMed

@dataclass
class DegValidation(AiResult):
    comparison: str = ""
    expected_genes: list[str] = field(default_factory=list)
    unexpected_genes: list[str] = field(default_factory=list)
    literature_links: list[GeneLitRef] = field(default_factory=list)
    validation_summary: str = ""
    discovery_highlights: list[str] = field(default_factory=list)

def run(
    adata,
    config: dict,
    study_context: dict,
    *,
    n_top_genes: int = 10,      # top N DEGs per comparison to validate
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> list[DegValidation] | None:
    ...
```

### What it reads from adata
- `adata.uns['rank_genes_groups']` — DEG results per cluster/comparison
- `adata.obs['ai_cell_type']` if available (from B1), else falls back to
  cluster ID as cell type label
- If rank_genes_groups missing → raise informative error

### Two-part analysis per comparison
1. **Validation (LLM)**: are top DEGs consistent with known biology?
   Returns expected_genes, unexpected_genes, validation_summary,
   discovery_highlights.
2. **Literature linking (PubMed RAG)**: for each top gene, retrieve
   PMID + title via PubMed search. One PubMed query per gene.
   Store as GeneLitRef. PMID + title only — never abstract text.

### PubMed query pattern per gene
```
"{gene} {tissue} {disease_context} single-cell"
```

### Copyright constraint (carry into every session)
Only PMIDs + titles stored and reported. Never abstract text, never
quotes from papers.

### Usage pattern (identical to all feature modules)
```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config, study_context, *, n_top_genes, log_dir, runtime_ai):
    try:
        check_ai_enabled(config, module="deg_validator", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None
```

## Step 5 — Tests (all without a real API key)
Mock pattern:
  - patch("ai.deg_validator.call_llm") — returns str only
  - patch("ai.deg_validator._query_pubmed") — returns list[GeneLitRef]
    (mock PubMed so tests never hit the network)

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns list[DegValidation] when mock LLM returns valid JSON
  - expected_genes / unexpected_genes populated correctly
  - literature_links populated from mock PubMed results
  - Missing rank_genes_groups raises informative error
  - Empty DEG list returns empty result gracefully
  - PMID deduplication across genes (same PMID for two genes → stored once)
  - AiResult base fields (timestamp, model, provider, skill_name) populated
  - Audit log written to log_dir after successful run
  - No abstract text stored in any GeneLitRef field

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ✅ DONE — A1: Pipeline advisor (13 tests passing)
Session 2  ✅ DONE — A2: Clustering advisor (22 tests passing)
Session 3  ✅ DONE — B1: Cluster annotator (23 tests passing)
Session 4  ← TODAY — B2: DEG validator + literature linker
Session 5  — B3: Coherence reviewer (build analysis_summary.json here)
Session 6  — A3: Downstream analysis suggester
Session 7  — C1: Narrative generator
Session 8  — C2: Full report + PowerPoint
Session 9  — Milestone validation: groundedness test + end-to-end GSE166635
```
Full details: .dev_memory/PHASE3_PLAN.md

## AI Layer Design Principles (carry into every session)
- Skills are the prompts — no raw strings in Python, ever
- AI is optional at three levels: global / per-module / per-step runtime flag
- A disabled module always returns None silently — never errors, never logs
- Every LLM call is audit-logged — no exceptions
- Graceful degradation — a failed AI module returns None, never crashes
- Context scales quality — more study_context = better AI output
- analysis_summary.json is load-bearing — design carefully in Session 5
- Copyright enforced at report layer — PMIDs + titles only, no abstract text

## Confirmed API Signatures (do not guess — use these exactly)

### call_llm (ai/_llm_client.py)
```python
call_llm(
    skill_name: str,
    inputs: dict,
    config: dict,
    *,
    log_dir: str,
    module: str | None,
    runtime_ai: bool,
) -> str                   # raw response string ONLY — no tuple
```

### write_audit_record (ai/_audit_log.py)
```python
write_audit_record(
    *,
    log_dir: str | Path,
    module: str,
    skill_version: str,
    model: str,
    provider: str,
    input_summary: dict,
    token_usage: dict | None,
    raw_response: str,
    parsed_output: dict | None,
    parse_success: bool,
) -> None
```

### check_ai_enabled (ai/_config_gate.py)
```python
check_ai_enabled(config: dict, module: str, runtime_ai: bool = True) -> None
# Raises AiDisabledError if: global off, module off, or runtime_ai=False
```

### AiResult base fields (ai/_base.py)
```python
timestamp: str      # ISO-8601 UTC, set automatically
model: str
provider: str
skill_name: str     # ← NOTE: skill_name, NOT skill
skill_version: str
reasoning: str
```

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-3 Modules (done, do not modify)
- ai/pipeline_advisor.py            ← 13 tests passing
- ai/skills/pipeline_advisor.yaml
- ai/clustering_advisor.py          ← 22 tests passing
- ai/skills/clustering_advisor.yaml
- ai/cluster_annotator.py           ← 23 tests passing
- ai/skills/cluster_annotator.yaml

## Bugs Fixed This Session (carry as warnings)
- AiResult base field is skill_name (not skill) — confirmed from _base.py
- Always upload _base.py (and any infrastructure file) at session start
  if there is any doubt about field names — do not guess from memory

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- Use `2>&1 | tee logs/<dataset>.log` to capture all output
- run_pipeline.py must be at repo root
- obs['cell_type_vote'] is the consensus column (not obs['cell_type'])
- rpy2 NOT used — do not attempt R integration
- gseapy.enrichr returns NO Overlap column in v1.1.3
- harmonypy ho.Z_corr is (n_cells, n_pcs) — do NOT transpose
- obs[batch_key] must be cast .astype(str) before pandas operations
- pseudobulk requires layers['counts'] (raw integers) — NOT layers['logcounts']
- GSE166635 pseudobulk is DISABLED — only 2 samples
- test_phase0_structure.py: TestCLI asserts cfg.provider == "ollama" (not "claude")
- YAML skill files: output_schema must always use block scalar (|)
  Never use bare type annotations: list[{...}], str | null, list[str]
  These are illegal YAML and will cause ScannerError at load time.
- biochatter==0.14.2 pinned — do NOT upgrade without re-running all AI tests
- call_llm returns str only — never unpack as tuple
- write_audit_record uses token_usage=dict (or None), not separate kwargs
- Provider and model always read from config — never hardcoded anywhere
- Default provider: "ollama", default model: "llama3"

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2
