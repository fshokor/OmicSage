"""
ai/clustering_advisor.py
------------------------
Phase 3 — Session 2 — A2: Clustering Advisor

Suggests Leiden clustering resolution for scRNA-seq data.
First Phase 3 module that uses PubMed RAG.

Public API
----------
run(adata, config, study_context, resolution_sweep_results, *, log_dir, runtime_ai)
    → ClusteringAdvice | None

Returns None silently when:
  - config.ai.features is False
  - config.ai.modules.clustering_advisor is False
  - runtime_ai is False
  - LLM call fails or returns unparseable JSON

Never raises. Never modifies adata.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._audit_log import write_audit_record
from ai._llm_client import call_llm
from ai._skill_loader import load_skill

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------

@dataclass
class PubMedRef:
    pmid: str
    title: str
    resolution_used: float | None
    tissue: str
    disease: str | None


@dataclass
class ClusteringAdvice(AiResult):
    suggested_resolution: float = 0.5
    resolution_range: tuple[float, float] = (0.3, 0.8)
    expected_n_clusters: int = 0
    literature_context: list[PubMedRef] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Rule-based helpers (no LLM — deterministic)
# ---------------------------------------------------------------------------

def _best_from_sweep(
    resolution_sweep_results: dict[float, float],
) -> tuple[float, float]:
    """Return (best_resolution, best_silhouette) from sweep dict.

    Falls back to (0.5, 0.0) if sweep is empty.
    """
    if not resolution_sweep_results:
        return 0.5, 0.0
    best_res = max(resolution_sweep_results, key=resolution_sweep_results.get)
    return best_res, resolution_sweep_results[best_res]


def _resolution_bias(n_cells: int, base_resolution: float) -> float:
    """Apply cell-count bias to the base resolution suggestion.

    < 1 000 cells  → bias toward lower resolution (fewer clusters)
    > 50 000 cells → bias toward higher resolution
    """
    if n_cells < 1_000:
        return max(0.1, base_resolution - 0.1)
    if n_cells > 50_000:
        return min(2.0, base_resolution + 0.2)
    return base_resolution


def _sweep_summary_str(resolution_sweep_results: dict[float, float]) -> str:
    """Human-readable summary for the skill prompt."""
    if not resolution_sweep_results:
        return "No sweep results available."
    parts = [
        f"{res}: {sil:.4f}"
        for res, sil in sorted(resolution_sweep_results.items())
    ]
    return ", ".join(parts)


def _sweep_range(resolution_sweep_results: dict[float, float]) -> tuple[float, float]:
    """Return (min_resolution, max_resolution) from sweep keys."""
    if not resolution_sweep_results:
        return (0.3, 0.8)
    keys = sorted(resolution_sweep_results.keys())
    return (keys[0], keys[-1])


# ---------------------------------------------------------------------------
# PubMed RAG
# ---------------------------------------------------------------------------

def _pubmed_query_string(tissue: str, disease_context: str | None) -> str:
    """Build PubMed query string per PHASE3_PLAN spec."""
    parts = ["Leiden clustering resolution", tissue]
    if disease_context:
        parts.append(disease_context)
    parts.append("single-cell RNA-seq")
    return " ".join(parts)


def _fetch_pubmed_refs(
    tissue: str,
    disease_context: str | None,
    config: dict,
) -> list[dict]:
    """Query PubMed via BioChatter retrieval.

    Returns a list of dicts with keys: pmid, title, resolution_used, tissue, disease.
    Returns empty list on any failure — never raises.

    Copyright rule: PMIDs + titles only. Abstract text is never stored.
    """
    query = _pubmed_query_string(tissue, disease_context)
    try:
        ai_cfg = config.get("ai", {})
        provider = ai_cfg.get("provider", "ollama")
        model = ai_cfg.get("model", "llama3")

        from ai._llm_client import _build_conversation
        conv = _build_conversation(provider=provider, model=model, config=ai_cfg)

        if hasattr(conv, "get_pubmed_results"):
            raw_results = conv.get_pubmed_results(query, num_results=5)
        else:
            logger.warning(
                "clustering_advisor: BioChatter conversation has no "
                "get_pubmed_results() — skipping PubMed RAG."
            )
            return []

        refs = []
        for item in raw_results or []:
            pmid = str(item.get("pmid") or item.get("uid") or "")
            title = str(item.get("title") or "")
            if not pmid or not title:
                continue
            refs.append({
                "pmid": pmid,
                "title": title,
                "resolution_used": None,
                "tissue": tissue,
                "disease": disease_context,
            })
        return refs

    except Exception as exc:  # noqa: BLE001
        logger.warning("clustering_advisor: PubMed RAG failed — %s", exc)
        return []


def _refs_to_prompt_str(refs: list[dict]) -> str:
    """Format refs for skill prompt. PMIDs + titles only (copyright rule)."""
    if not refs:
        return ""
    parts = []
    for r in refs:
        resolution_str = (
            str(r["resolution_used"]) if r["resolution_used"] is not None else "not reported"
        )
        parts.append(
            f"{r['pmid']}|{r['title']}|{resolution_str}|{r['tissue']}|{r.get('disease') or 'N/A'}"
        )
    return "; ".join(parts)


# ---------------------------------------------------------------------------
# JSON parse helpers
# ---------------------------------------------------------------------------

def _parse_llm_response(raw: str, sweep_range: tuple[float, float]) -> ClusteringAdvice | None:
    """Parse LLM JSON response into ClusteringAdvice.

    Returns None on any parse failure — logs a warning but never raises.
    """
    import json

    try:
        cleaned = raw.strip()
        if cleaned.startswith("```"):
            cleaned = "\n".join(cleaned.split("\n")[1:])
        if cleaned.endswith("```"):
            cleaned = "\n".join(cleaned.split("\n")[:-1])
        data = json.loads(cleaned)
    except (json.JSONDecodeError, ValueError) as exc:
        logger.warning("clustering_advisor: JSON parse failed — %s", exc)
        return None

    try:
        suggested = float(data["suggested_resolution"])
        rng_raw = data["resolution_range"]
        resolution_range = (float(rng_raw[0]), float(rng_raw[1]))
        expected_n = int(data["expected_n_clusters"])
        reasoning = str(data.get("reasoning", ""))

        lit = []
        for item in data.get("literature_context", []):
            try:
                res_used = item.get("resolution_used")
                lit.append(
                    PubMedRef(
                        pmid=str(item["pmid"]),
                        title=str(item["title"]),
                        resolution_used=float(res_used) if res_used is not None else None,
                        tissue=str(item.get("tissue", "")),
                        disease=item.get("disease"),
                    )
                )
            except (KeyError, TypeError, ValueError) as e:
                logger.warning("clustering_advisor: skipping malformed ref — %s", e)

        # Sanity-check: resolution_range must be ordered
        if resolution_range[0] >= resolution_range[1]:
            lo = min(resolution_range)
            hi = max(resolution_range)
            resolution_range = (lo, hi) if lo != hi else (lo, lo + 0.2)

        return ClusteringAdvice(
            suggested_resolution=suggested,
            resolution_range=resolution_range,
            expected_n_clusters=expected_n,
            literature_context=lit,
            reasoning=reasoning,
        )

    except (KeyError, TypeError, ValueError) as exc:
        logger.warning("clustering_advisor: response field parse failed — %s", exc)
        return None


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------

def run(
    adata,
    config: dict,
    study_context: dict,
    resolution_sweep_results: dict[float, float],
    *,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> ClusteringAdvice | None:
    """Suggest Leiden clustering resolution using silhouette sweep + PubMed RAG.

    Parameters
    ----------
    adata
        AnnData after dimensionality reduction. Read-only — never modified.
    config
        Full pipeline config dict (must contain 'ai' section).
    study_context
        Loaded study_context.yaml as dict.
    resolution_sweep_results
        Dict mapping resolution (float) → silhouette score (float).
        Can be empty — a default resolution is used with a warning.
    log_dir
        Directory for audit log JSONL files.
    runtime_ai
        Set False via --no-ai flag at runtime to skip all AI calls.

    Returns
    -------
    ClusteringAdvice or None
        None if AI is disabled, sweep is unresolvable, or LLM call fails.
    """
    # --- Gate check (three-level: global / per-module / runtime flag) -------
    try:
        check_ai_enabled(config, module="clustering_advisor", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # --- Study context fields -----------------------------------------------
    tissue = study_context.get("dataset", {}).get("tissue", "unknown tissue")
    disease_context = (
        study_context.get("disease", {}).get("context") or "not applicable"
    )
    n_cells = adata.n_obs if adata is not None else 0
    n_hvg = (
        int(adata.var["highly_variable"].sum())
        if adata is not None and "highly_variable" in adata.var.columns
        else 0
    )

    # --- Rule-based pre-checks ----------------------------------------------
    if not resolution_sweep_results:
        warnings.warn(
            "clustering_advisor: resolution_sweep_results is empty — "
            "using default resolution 0.5.",
            UserWarning,
            stacklevel=2,
        )

    best_res, best_sil = _best_from_sweep(resolution_sweep_results)
    biased_res = _resolution_bias(n_cells, best_res)
    sweep_summary = _sweep_summary_str(resolution_sweep_results)
    sweep_rng = _sweep_range(resolution_sweep_results)

    # --- PubMed RAG ---------------------------------------------------------
    pubmed_refs = _fetch_pubmed_refs(tissue, disease_context, config)
    literature_str = _refs_to_prompt_str(pubmed_refs)

    # --- Build skill inputs -------------------------------------------------
    skill_inputs = {
        "tissue": tissue,
        "disease_context": disease_context,
        "n_cells": n_cells,
        "n_highly_variable_genes": n_hvg,
        "best_resolution": biased_res,
        "best_silhouette": round(best_sil, 4),
        "resolution_sweep_summary": sweep_summary,
        "literature_refs": literature_str if literature_str else "No references found.",
    }

    # --- LLM call -----------------------------------------------------------
    # call_llm signature: (skill_name, inputs, config, *, log_dir, module)
    # Returns: raw response string only — provider/model come from config
    ai_cfg = config.get("ai", {})
    provider = ai_cfg.get("provider", "ollama")
    model = ai_cfg.get("model", "llama3")

    raw_response = call_llm(
        skill_name="clustering_advisor",
        inputs=skill_inputs,
        config=config,
        log_dir=log_dir,
        module="clustering_advisor",
    )

    # --- Parse --------------------------------------------------------------
    advice = _parse_llm_response(raw_response, sweep_rng)

    # --- Audit log ----------------------------------------------------------
    # call_llm already writes a partial audit record; we write a second record
    # with the parsed output so the final parse result is captured.
    def _advice_to_dict(a: ClusteringAdvice) -> dict:
        return {
            "suggested_resolution": a.suggested_resolution,
            "resolution_range": list(a.resolution_range),
            "expected_n_clusters": a.expected_n_clusters,
            "literature_context": [r.__dict__ for r in a.literature_context],
            "reasoning": a.reasoning,
        }

    write_audit_record(
        log_dir=log_dir,
        module="clustering_advisor",
        skill_version="1.0",
        model=model,
        provider=provider,
        input_summary={
            "tissue": tissue,
            "disease_context": disease_context,
            "n_cells": n_cells,
            "n_hvg": n_hvg,
            "best_resolution": biased_res,
            "n_pubmed_refs": len(pubmed_refs),
        },
        token_usage=None,
        raw_response=raw_response,
        parsed_output=_advice_to_dict(advice) if advice else None,
        parse_success=advice is not None,
    )

    if advice is None:
        return None

    # Stamp base fields from config (what was actually used)
    import datetime
    advice.timestamp = datetime.datetime.utcnow().isoformat() + "Z"
    advice.model = model
    advice.provider = provider

    return advice
