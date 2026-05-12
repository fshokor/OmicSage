## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Full Phase 3 plan designed and saved to
                      .dev_memory/PHASE3_PLAN.md. Three foundation
                      files built and tested:
                        - config/study_context_template.yaml
                        - ai/skills/cluster_annotator.yaml  (reference pattern)
                        - ai/_skill_loader.py               (tested, all checks pass)
File last worked on: ai/_skill_loader.py

## Today's Goal
Phase 3 — Session 0: Shared AI infrastructure.

Build the four foundation modules that every AI feature will import.
Nothing else. No feature modules today.

Files to create:
  - ai/_llm_client.py      ← BioChatter wrapper, routes on provider
  - ai/_audit_log.py       ← write_audit_record() → logs/llm/<module>.jsonl
  - ai/_config_gate.py     ← raises AiDisabledError if ai_features: false
  - ai/_base.py            ← AiResult base dataclass
  - tests/test_ai_infrastructure.py

## Step 0 — Read the Phase 3 plan before writing any code
```
.dev_memory/PHASE3_PLAN.md
```
The full plan is there. Every session references it.

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~231 passed, 1-2 skipped
```

## Step 2 — Verify Phase 2 milestone still works
```bash
python run_pipeline.py --config config/runs/GSE166635.yaml --step qc
# Should cache and auto-generate reports/GSE166635/00_combined_report.html
```

## Step 3 — Install and verify biochatter
```bash
conda activate omicsage
pip install biochatter==0.14.2   # pin this exact version
python -c "
from biochatter.llm_connect import OllamaConversation
conv = OllamaConversation(
    base_url='http://localhost:11434',
    prompts={}, model_name='llama3', correct=False,
)
conv.append_system_message('You are a helpful assistant.')
response, _, _ = conv.query('Say hello in one sentence.')
print(response)
"
```

BioChatter API verified on v0.14.2 — use these exact method names in _llm_client.py:
  - append_system_message(text)   ← NOT set_system_message (renamed in v0.14)
  - query(text)                   → returns (response, token_usage, correction)
  - set_api_key(api_key=...)      ← Claude and OpenAI only, not Ollama
  - base_url is required for OllamaConversation ← added in recent version

## Step 4 — Copy foundation files already built this planning session
These files were designed and tested — copy them into the repo before
writing any new code:

```bash
# Skill loader (already tested — all checks pass)
cp <delivered_files>/ai/_skill_loader.py ai/_skill_loader.py

# Reference skill file (establishes pattern for all other skills)
mkdir -p ai/skills
cp <delivered_files>/ai/skills/cluster_annotator.yaml ai/skills/cluster_annotator.yaml

# Study context template
cp <delivered_files>/config/study_context_template.yaml config/study_context_template.yaml

# Populate study context for each existing dataset
cp config/study_context_template.yaml config/runs/GSE166635/study_context.yaml
cp config/study_context_template.yaml config/runs/GSE194122/study_context.yaml
# Then fill in tissue / disease / experiment fields for each dataset
```

## Step 5 — Build the four infrastructure modules

### ai/_config_gate.py
Three-level optionality — all checked by a single function:
  Level 1 — global:     ai_features: false → disabled
  Level 2 — per-module: modules: { cluster_annotator: false } → disabled
  Level 3 — runtime:    --ai flag absent on run_pipeline.py → disabled

```python
# Usage pattern in every feature module — identical every time:
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config):
    try:
        check_ai_enabled(config, module="cluster_annotator")
    except AiDisabledError:
        return None
    # ... rest of the module
```

A module that is off returns None silently. Never raises. Never logs a warning.
Real errors (network, parse failure, bad input) propagate normally — never caught here.

### ai/_base.py
AiResult base dataclass. All feature dataclasses inherit from it.
Fields: timestamp (ISO-8601 UTC), model, provider, skill_name,
        skill_version, reasoning.

### ai/_audit_log.py
write_audit_record() appends one JSONL line per LLM call.
Record format (from PHASE3_PLAN.md):
  timestamp, module, skill_version, model, provider,
  input_summary, raw_response, parsed_output, parse_success

### ai/_llm_client.py
Thin wrapper around BioChatter. Routes on config.ai.provider.
Single public function: call_llm(skill_name, inputs, config) → str
  - Loads skill via _skill_loader.load_skill()
  - Routes to AnthropicConversation | OllamaConversation | GptConversation
  - Writes audit log via _audit_log.write_audit_record()
  - Returns raw LLM response string
  - Feature modules handle JSON parsing themselves

Config section expected:
```yaml
ai:
  features: true
  provider: claude        # claude | ollama | openai
  model: claude-sonnet-4-20250514
```

## Infrastructure Tests Required (all without a real API key)
- Config gate raises AiDisabledError when ai_features: false
- Config gate raises AiDisabledError when module explicitly set to false
- Config gate passes when module key is missing (defaults to enabled)
- Config gate raises when --ai flag absent at runtime (runtime_ai=False)
- LLM client routes to AnthropicConversation when provider: claude
- LLM client routes to OllamaConversation when provider: ollama
- LLM client routes to GptConversation when provider: openai
- Unknown provider raises ValueError listing valid options
- Audit log creates file, writes correct JSONL, appends on second call
- Audit log creates log_dir if it does not exist
- AiResult base dataclass instantiates correctly
- Skill loader integrates with LLM client (mock provider)

Hard stop: do not start ai/pipeline_advisor.py until all infrastructure
tests pass.

## Phase 3 Build Order (full reference)
```
Session 0  ← TODAY — Infrastructure: _llm_client, _audit_log,
                      _skill_loader (done), _config_gate, _base
Session 1  — A1: Pipeline advisor
Session 2  — A2: Clustering advisor (first PubMed RAG use)
Session 3  — B1: Cluster annotator
Session 4  — B2: DEG validator + literature linker
Session 5  — B3: Coherence reviewer (build analysis_summary.json here)
Session 6  — A3: Downstream analysis suggester
Session 7  — C1: Narrative generator
Session 8  — C2: Full report + PowerPoint
Session 9  — Milestone validation: groundedness test + end-to-end GSE166635
```
Full details for every session: .dev_memory/PHASE3_PLAN.md

## AI Layer Design Principles (carry into every session)
- Skills are the prompts — no raw strings in Python, ever
- AI is optional at three levels: global / per-module / per-step runtime flag
- A disabled module always returns None silently — never errors, never logs
- Every LLM call is audit-logged — no exceptions
- Graceful degradation — a failed AI module returns None, never crashes
- Context scales quality — more study_context = better AI output
- analysis_summary.json is load-bearing — design carefully in Session 5
- Copyright enforced at report layer — PMIDs + titles only, no abstract text

## New Files Created This Planning Session
- ai/_skill_loader.py                       ← tested, ready to copy in
- ai/skills/cluster_annotator.yaml          ← reference skill pattern
- config/study_context_template.yaml        ← copy into each run dir
- .dev_memory/PHASE3_PLAN.md                ← full Phase 3 plan

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

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2
New dependency this phase: biochatter (pip install biochatter)
