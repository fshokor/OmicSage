"""
ai/downstream_suggester.py
Phase 3 — Session 6: A3 Downstream Analysis Suggester

Suggests follow-up analysis steps based on what was found in the completed
single-cell analysis. Rule-based triggers fire first (no LLM needed); the LLM
then adds context-aware rationale and additional suggestions.

Public API
----------
run(adata, config, study_context, coherence_review=None, *,
    output_path=None, log_dir="logs/llm", runtime_ai=True)
    -> DownstreamAdvice | None
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class DownstreamSuggestion:
    step_name: str = ""
    rationale: str = ""
    expected_output: str = ""
    relevant_tool: str = ""


@dataclass
class DownstreamAdvice(AiResult):
    suggestions: list[DownstreamSuggestion] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Rule-based trigger helpers (no LLM — always evaluated)
# ---------------------------------------------------------------------------

# Cell type keywords used to classify clusters as progenitor / mature / immune / non-immune
_PROGENITOR_KEYWORDS = {
    "progenitor", "stem", "precursor", "blast", "progenitors",
    "hsc", "lsc", "cmp", "gmp", "mep",
}
_MATURE_KEYWORDS = {
    "mature", "differentiated", "effector", "memory", "naive",
    "macrophage", "monocyte", "neutrophil", "nk", "b cell", "t cell",
    "hepatocyte", "fibroblast", "endothelial", "epithelial",
}
_IMMUNE_KEYWORDS = {
    "t cell", "b cell", "nk", "macrophage", "monocyte", "dendritic",
    "neutrophil", "mast", "plasma", "lymphocyte", "cd8", "cd4",
    "treg", "nkt", "myeloid", "lymphoid",
}
_NON_IMMUNE_KEYWORDS = {
    "hepatocyte", "fibroblast", "endothelial", "epithelial", "stellate",
    "cholangiocyte", "tumour", "cancer", "malignant", "stromal",
    "smooth muscle", "pericyte", "astrocyte", "neuron", "cardiomyocyte",
}


def _cell_type_labels(adata) -> list[str]:
    """Return unique cell type labels from obs, lowercased."""
    for col in ("ai_cell_type", "cell_type_vote", "cell_type"):
        if col in adata.obs.columns:
            return [str(v).lower() for v in adata.obs[col].unique()]
    return []


def _matches_any(label: str, keywords: set[str]) -> bool:
    return any(kw in label for kw in keywords)


def _has_progenitor_and_mature(labels: list[str]) -> bool:
    return any(_matches_any(l, _PROGENITOR_KEYWORDS) for l in labels) and \
           any(_matches_any(l, _MATURE_KEYWORDS) for l in labels)


def _has_immune_and_non_immune(labels: list[str]) -> bool:
    return any(_matches_any(l, _IMMUNE_KEYWORDS) for l in labels) and \
           any(_matches_any(l, _NON_IMMUNE_KEYWORDS) for l in labels)


def _has_clinical_metadata(adata) -> bool:
    clinical_cols = {
        "survival", "os", "pfs", "stage", "grade", "response",
        "treatment", "recurrence", "event", "time_to_event",
    }
    obs_cols = {c.lower() for c in adata.obs.columns}
    return bool(obs_cols & clinical_cols)


def _rule_based_suggestions(
    adata,
    study_context: dict,
    coherence_review: Any | None,
) -> list[DownstreamSuggestion]:
    """
    Evaluate rule-based triggers and return pre-formed suggestions.
    These always run — no LLM required.
    """
    suggestions: list[DownstreamSuggestion] = []
    labels = _cell_type_labels(adata)

    # Rule 1: progenitor + mature → trajectory
    if _has_progenitor_and_mature(labels):
        suggestions.append(DownstreamSuggestion(
            step_name="Trajectory analysis",
            rationale=(
                "Both progenitor and mature cell types are present in the dataset, "
                "suggesting a differentiation continuum that trajectory analysis can resolve."
            ),
            expected_output=(
                "Pseudotime ordering of cells along differentiation axes, "
                "with branch points identifying fate decisions."
            ),
            relevant_tool="Slingshot / PAGA (scVelo for RNA velocity if spliced counts available)",
        ))

    # Rule 2: immune + non-immune → cell-cell communication
    if _has_immune_and_non_immune(labels):
        suggestions.append(DownstreamSuggestion(
            step_name="Cell-cell communication analysis",
            rationale=(
                "Immune and non-immune cell populations co-exist in the dataset. "
                "Ligand-receptor analysis will reveal signalling interactions "
                "between the microenvironment compartments."
            ),
            expected_output=(
                "Ranked ligand-receptor pairs between cell type pairs, "
                "with pathway-level summaries of dominant signalling axes."
            ),
            relevant_tool="CellChat / LIANA",
        ))

    # Rule 3: clinical metadata → survival analysis
    if _has_clinical_metadata(adata):
        suggestions.append(DownstreamSuggestion(
            step_name="Survival analysis",
            rationale=(
                "Clinical outcome metadata detected in obs. Cell type abundances "
                "or gene signatures can be correlated with patient survival."
            ),
            expected_output=(
                "Kaplan-Meier curves stratified by cell type abundance or "
                "gene signature score, with log-rank p-values."
            ),
            relevant_tool="lifelines / survminer (R)",
        ))

    # Rule 4: coherence review flags sub-clustering candidates
    if coherence_review is not None:
        candidates = getattr(coherence_review, "sub_clustering_candidates", [])
        if candidates:
            candidate_str = ", ".join(str(c) for c in candidates)
            suggestions.append(DownstreamSuggestion(
                step_name="Sub-clustering of heterogeneous populations",
                rationale=(
                    f"Coherence review identified cluster(s) {candidate_str} as "
                    "containing mixed cell type markers, suggesting unresolved "
                    "heterogeneity that finer resolution clustering can reveal."
                ),
                expected_output=(
                    "Sub-clusters within the flagged populations, with marker "
                    "genes distinguishing each sub-population."
                ),
                relevant_tool="Leiden clustering at higher resolution (scanpy.tl.leiden)",
            ))

    # Rule 5: multiple conditions → pseudobulk DEG (if not already done)
    n_conditions = study_context.get("experiment", {}).get("n_conditions", 1)
    if n_conditions and int(n_conditions) > 1:
        # Check whether pseudobulk was already run
        already_pseudobulk = False
        if hasattr(adata, "uns") and "pseudobulk_deg" in adata.uns:
            already_pseudobulk = True
        if not already_pseudobulk:
            suggestions.append(DownstreamSuggestion(
                step_name="Pseudobulk differential expression",
                rationale=(
                    f"The dataset has {n_conditions} conditions. Pseudobulk DEG "
                    "aggregates counts per donor before testing, providing better "
                    "control of donor-level variance than single-cell Wilcoxon tests."
                ),
                expected_output=(
                    "Per-condition DEG tables with effect sizes and adjusted p-values, "
                    "suitable for downstream pathway enrichment."
                ),
                relevant_tool="PyDESeq2 / DESeq2 (R)",
            ))

    return suggestions


# ---------------------------------------------------------------------------
# NEXT_STEPS.md writer
# ---------------------------------------------------------------------------

def _write_next_steps_md(
    suggestions: list[DownstreamSuggestion],
    reasoning: str,
    output_path: str | Path,
) -> None:
    lines = [
        "# OmicSage — Suggested Next Steps\n",
        f"*Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}*\n",
        "",
        "## Overview\n",
        reasoning,
        "",
    ]
    for i, s in enumerate(suggestions, 1):
        lines += [
            f"## {i}. {s.step_name}\n",
            f"**Rationale**: {s.rationale}\n",
            f"**Expected output**: {s.expected_output}\n",
            f"**Relevant tool**: {s.relevant_tool}\n",
            "",
        ]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")
    logger.info("NEXT_STEPS.md written to %s", output_path)


# ---------------------------------------------------------------------------
# JSON parsing helper
# ---------------------------------------------------------------------------

def _parse_response(raw: str) -> tuple[list[DownstreamSuggestion], str]:
    """
    Parse LLM JSON response into (suggestions, reasoning).
    Returns ([], "") on any parse failure — never raises.
    """
    try:
        # Strip markdown fences if present
        text = raw.strip()
        if text.startswith("```"):
            lines = text.splitlines()
            text = "\n".join(lines[1:-1] if lines[-1].strip() == "```" else lines[1:])
        data = json.loads(text)
        suggestions = []
        for item in data.get("suggestions", []):
            suggestions.append(DownstreamSuggestion(
                step_name=item.get("step_name", ""),
                rationale=item.get("rationale", ""),
                expected_output=item.get("expected_output", ""),
                relevant_tool=item.get("relevant_tool", ""),
            ))
        reasoning = data.get("reasoning", "")
        return suggestions, reasoning
    except Exception as exc:
        logger.warning("downstream_suggester: failed to parse LLM response — %s", exc)
        return [], ""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(
    adata,
    config: dict,
    study_context: dict,
    coherence_review=None,
    *,
    output_path: str | None = None,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> DownstreamAdvice | None:
    """
    Suggest downstream analyses given completed single-cell analysis findings.

    Parameters
    ----------
    adata:
        Annotated data matrix (AnnData). Used for rule-based trigger evaluation.
    config:
        Pipeline config dict. Must contain ai section.
    study_context:
        Loaded study_context.yaml as dict.
    coherence_review:
        CoherenceReview dataclass from B3, or None.
    output_path:
        If set, NEXT_STEPS.md is written to this path.
    log_dir:
        Directory for LLM audit logs.
    runtime_ai:
        False when --ai flag absent from CLI. Causes immediate return None.

    Returns
    -------
    DownstreamAdvice or None (if AI disabled at any level).
    """
    try:
        check_ai_enabled(config, module="downstream_suggester", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # --- Rule-based triggers (always evaluated before LLM) ---
    rule_suggestions = _rule_based_suggestions(adata, study_context, coherence_review)

    # --- Build LLM inputs ---
    tissue = study_context.get("dataset", {}).get("tissue", "unknown")
    disease_context = study_context.get("disease", {}).get("context", None)
    biological_question = study_context.get("biological_question", "")

    # Serialise analysis_summary if available from adata.uns
    analysis_summary_dict = {}
    if hasattr(adata, "uns") and "omicsage_analysis_summary" in adata.uns:
        analysis_summary_dict = adata.uns["omicsage_analysis_summary"]
    analysis_summary_str = json.dumps(analysis_summary_dict, indent=2)

    # Serialise coherence flags
    coherence_flags = []
    if coherence_review is not None:
        raw_flags = getattr(coherence_review, "flags", [])
        for f in raw_flags:
            if hasattr(f, "__dict__"):
                coherence_flags.append(f.__dict__)
            elif isinstance(f, dict):
                coherence_flags.append(f)
    coherence_flags_str = json.dumps(coherence_flags, indent=2)

    inputs = {
        "tissue": tissue,
        "disease_context": str(disease_context) if disease_context else "none",
        "biological_question": biological_question or "Not specified",
        "analysis_summary": analysis_summary_str,
        "coherence_flags": coherence_flags_str,
    }

    # --- LLM call ---
    raw_response = call_llm(
        skill_name="downstream_suggester",
        inputs=inputs,
        config=config,
        log_dir=log_dir,
        module="downstream_suggester",
        runtime_ai=runtime_ai,
    )

    llm_suggestions, reasoning = _parse_response(raw_response)

    # Merge: rule-based first, then LLM additions (deduplicate by step_name)
    existing_names = {s.step_name.lower() for s in rule_suggestions}
    for s in llm_suggestions:
        if s.step_name.lower() not in existing_names:
            rule_suggestions.append(s)
            existing_names.add(s.step_name.lower())

    all_suggestions = rule_suggestions

    # Fallback reasoning if LLM parse failed
    if not reasoning:
        reasoning = (
            "Rule-based analysis identified potential downstream steps based on "
            "detected cell types, experimental design, and coherence review flags."
        )

    # --- Write NEXT_STEPS.md ---
    if output_path is not None:
        _write_next_steps_md(all_suggestions, reasoning, output_path)

    # --- Build result ---
    ai_cfg = config.get("ai", {})
    return DownstreamAdvice(
        timestamp=datetime.now(timezone.utc).isoformat(),
        model=ai_cfg.get("model", "llama3"),
        provider=ai_cfg.get("provider", "ollama"),
        skill_name="downstream_suggester",
        skill_version="1.0",
        reasoning=reasoning,
        suggestions=all_suggestions,
    )
