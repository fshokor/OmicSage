"""
ai/deg_validator.py
-------------------
Phase 3 — Session 4: B2 — DEG Validator + Literature Linker

Two-part analysis per DEG comparison:
  1. LLM validation — are top DEGs consistent with known biology?
  2. PubMed RAG — retrieve PMID + title per gene (no abstract text ever)

Public API
----------
  run(adata, config, study_context, *, n_top_genes, log_dir, runtime_ai)
      -> list[DegValidation] | None

Returned None means AI is disabled (global, module, or runtime flag).
Raises ValueError if adata.uns['rank_genes_groups'] is missing.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm
from ai._skill_loader import load_skill
from ai._audit_log import write_audit_record

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------

@dataclass
class GeneLitRef:
    """One PubMed reference linked to a gene. PMID + title only — never abstract."""
    gene: str = ""
    pmid: str = ""
    title: str = ""
    context: str = ""   # one-sentence context derived from title, not from abstract


@dataclass
class DegValidation(AiResult):
    """Validation result for one DEG comparison."""
    comparison: str = ""
    expected_genes: list[str] = field(default_factory=list)
    unexpected_genes: list[str] = field(default_factory=list)
    literature_links: list[GeneLitRef] = field(default_factory=list)
    validation_summary: str = ""
    discovery_highlights: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# PubMed helper (mockable in tests)
# ---------------------------------------------------------------------------

def _query_pubmed(
    gene: str,
    tissue: str,
    disease_context: str | None,
) -> list[GeneLitRef]:
    """
    Query PubMed for one gene and return up to 3 GeneLitRef results.

    Query pattern: "{gene} {tissue} {disease_context} single-cell"

    Returns PMIDs + titles only. Never stores abstract text.
    Falls back to empty list on any network or parse error.
    """
    try:
        from Bio import Entrez  # biopython — available in omicsage env
    except ImportError:
        logger.warning("biopython not installed — PubMed queries skipped")
        return []

    disease_part = disease_context if disease_context else ""
    query = f"{gene} {tissue} {disease_part} single-cell".strip()

    try:
        Entrez.email = "omicsage@example.com"
        # Search phase
        handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
        search_record = Entrez.read(handle)
        handle.close()

        pmids = search_record.get("IdList", [])
        if not pmids:
            return []

        # Fetch titles only — no abstract, no full text
        handle = Entrez.efetch(
            db="pubmed",
            id=",".join(pmids),
            rettype="xml",
            retmode="xml",
        )
        fetch_records = Entrez.read(handle)
        handle.close()

        results: list[GeneLitRef] = []
        for record in fetch_records.get("PubmedArticle", []):
            try:
                pmid = str(
                    record["MedlineCitation"]["PMID"]
                )
                article = record["MedlineCitation"]["Article"]
                title = str(article.get("ArticleTitle", ""))
                # One-sentence context derived from title alone — never from abstract
                context = f"{gene} mentioned in: {title}"
                results.append(
                    GeneLitRef(gene=gene, pmid=pmid, title=title, context=context)
                )
            except (KeyError, TypeError):
                continue

        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("PubMed query failed for gene %s: %s", gene, exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_degs(
    adata: Any,
    n_top_genes: int,
) -> dict[str, list[dict]]:
    """
    Extract top N DEGs per group from adata.uns['rank_genes_groups'].

    Returns dict: group_name -> list of dicts with keys
      gene, log2fc, padj, cell_type
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "adata.uns['rank_genes_groups'] is missing. "
            "Run sc.tl.rank_genes_groups() before calling the DEG validator."
        )

    rgg = adata.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)
    result: dict[str, list[dict]] = {}

    for group in groups:
        names = [rgg["names"][group][i] for i in range(min(n_top_genes, len(rgg["names"][group])))]
        # log2fc and pvals_adj may not always be present (depends on method)
        log2fcs = (
            [float(rgg["logfoldchanges"][group][i]) for i in range(len(names))]
            if "logfoldchanges" in rgg
            else [0.0] * len(names)
        )
        padjs = (
            [float(rgg["pvals_adj"][group][i]) for i in range(len(names))]
            if "pvals_adj" in rgg
            else [1.0] * len(names)
        )

        # Cell type per cell — use ai_cell_type if available, else group label
        cell_type: str | None = None
        if "ai_cell_type" in adata.obs.columns:
            # Most common ai_cell_type in this group
            mask = adata.obs["leiden"] == str(group) if "leiden" in adata.obs.columns else None
            if mask is not None and mask.any():
                cell_type = adata.obs.loc[mask, "ai_cell_type"].mode().iloc[0]

        degs = [
            {
                "gene": names[i],
                "log2fc": log2fcs[i],
                "padj": padjs[i],
                "cell_type": cell_type,
            }
            for i in range(len(names))
        ]
        result[group] = degs

    return result


def _parse_llm_response(raw: str) -> dict:
    """
    Parse LLM JSON response. Strips markdown fences if present.
    Returns dict with expected keys or raises ValueError.
    """
    clean = raw.strip()
    if not clean:
        raise ValueError("LLM returned an empty response")

    if clean.startswith("```"):
        lines = clean.splitlines()
        # Drop first and last fence lines
        clean = "\n".join(lines[1:-1] if lines[-1].strip() == "```" else lines[1:])

    parsed = json.loads(clean)

    required_keys = {
        "expected_genes",
        "unexpected_genes",
        "validation_summary",
        "discovery_highlights",
    }
    missing = required_keys - set(parsed.keys())
    if missing:
        raise ValueError(f"LLM response missing keys: {missing}")

    return parsed


def _deduplicate_lit_refs(refs: list[GeneLitRef]) -> list[GeneLitRef]:
    """
    Deduplicate GeneLitRef list by PMID. When the same PMID appears for
    multiple genes, keep the first occurrence and append additional gene
    names to the context field.
    """
    seen: dict[str, GeneLitRef] = {}
    for ref in refs:
        if ref.pmid not in seen:
            seen[ref.pmid] = ref
        else:
            # Extend context to mention the additional gene
            existing = seen[ref.pmid]
            if ref.gene not in existing.context:
                existing.context = f"{existing.context}; also relevant for {ref.gene}"
    return list(seen.values())


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(
    adata: Any,
    config: dict,
    study_context: dict,
    *,
    n_top_genes: int = 10,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> list[DegValidation] | None:
    """
    Validate DEGs against known biology and link to PubMed literature.

    Parameters
    ----------
    adata
        AnnData with adata.uns['rank_genes_groups'] populated.
    config
        Full pipeline config dict (must contain ai section).
    study_context
        Contents of study_context.yaml as a dict.
    n_top_genes
        Number of top DEGs per group to validate.
    log_dir
        Directory for audit log JSONL files.
    runtime_ai
        Set to False to disable AI for this call regardless of config.

    Returns
    -------
    list[DegValidation] or None
        None if AI is disabled at any level.
        Empty list if rank_genes_groups has no groups.
    """
    try:
        check_ai_enabled(config, module="deg_validator", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    tissue = study_context.get("dataset", {}).get("tissue", "unknown")
    disease_context = study_context.get("disease", {}).get("context", None)

    # Extract DEGs — raises ValueError if rank_genes_groups missing
    groups_degs = _extract_degs(adata, n_top_genes)

    if not groups_degs:
        return []

    # Load skill metadata for audit log
    skill_meta = load_skill("deg_validator", tissue="", disease_context="", comparison="", degs=[])
    skill_version = "1.0"  # read from skill file at load; hardcoded as fallback

    results: list[DegValidation] = []

    for group, degs in groups_degs.items():
        if not degs:
            # Empty DEG list for this group — return empty DegValidation
            results.append(
                DegValidation(
                    timestamp=datetime.now(timezone.utc).isoformat(),
                    model=config.get("ai", {}).get("model", "llama3"),
                    provider=config.get("ai", {}).get("provider", "ollama"),
                    skill_name="deg_validator",
                    skill_version=skill_version,
                    reasoning="No DEGs provided for this group.",
                    comparison=group,
                )
            )
            continue

        # ------------------------------------------------------------------
        # Part 1: LLM validation
        # ------------------------------------------------------------------
        inputs = {
            "tissue": tissue,
            "disease_context": disease_context if disease_context else "null",
            "comparison": group,
            "degs": json.dumps(degs, indent=2),
        }

        parse_success = False
        parsed: dict = {}
        last_exc: Exception | None = None

        for attempt in range(1, 3):   # up to 2 attempts
            raw_response = call_llm(
                skill_name="deg_validator",
                inputs=inputs,
                config=config,
                log_dir=log_dir,
                module="deg_validator",
                runtime_ai=runtime_ai,
            )
            try:
                parsed = _parse_llm_response(raw_response)
                parse_success = True
                break
            except (json.JSONDecodeError, ValueError) as exc:
                last_exc = exc
                if attempt < 2:
                    logger.warning(
                        "DEG validator: attempt %d failed for group %s: %s — retrying",
                        attempt,
                        group,
                        exc,
                    )

        if not parse_success:
            logger.warning(
                "DEG validator: failed to parse LLM response for group %s after 2 attempts: %s",
                group,
                last_exc,
            )

        # ------------------------------------------------------------------
        # Part 2: PubMed RAG — one query per gene
        # ------------------------------------------------------------------
        all_refs: list[GeneLitRef] = []
        gene_names = [d["gene"] for d in degs]

        for gene in gene_names:
            refs = _query_pubmed(gene, tissue, disease_context)
            all_refs.extend(refs)

        # Deduplicate by PMID
        deduped_refs = _deduplicate_lit_refs(all_refs)

        # ------------------------------------------------------------------
        # Audit log
        # ------------------------------------------------------------------
        model = config.get("ai", {}).get("model", "llama3")
        provider = config.get("ai", {}).get("provider", "ollama")

        write_audit_record(
            log_dir=log_dir,
            module="deg_validator",
            skill_version=skill_version,
            model=model,
            provider=provider,
            input_summary={
                "tissue": tissue,
                "disease_context": disease_context,
                "comparison": group,
                "n_degs": len(degs),
                "n_genes": len(gene_names),
            },
            token_usage=None,
            raw_response=raw_response,
            parsed_output=parsed if parse_success else None,
            parse_success=parse_success,
        )

        # ------------------------------------------------------------------
        # Assemble result
        # ------------------------------------------------------------------
        results.append(
            DegValidation(
                timestamp=datetime.now(timezone.utc).isoformat(),
                model=model,
                provider=provider,
                skill_name="deg_validator",
                skill_version=skill_version,
                reasoning=parsed.get("validation_summary", ""),
                comparison=group,
                expected_genes=parsed.get("expected_genes", []),
                unexpected_genes=parsed.get("unexpected_genes", []),
                literature_links=deduped_refs,
                validation_summary=parsed.get("validation_summary", ""),
                discovery_highlights=parsed.get("discovery_highlights", []),
            )
        )

    return results
