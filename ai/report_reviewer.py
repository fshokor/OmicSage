"""
ai/report_reviewer.py
---------------------
D1 — Report Reviewer

Reads the final HTML report, strips tags, calls the LLM to review scientific
quality, and writes a structured report_review.md to the run output directory.

Public API
----------
run(
    html_report_path: str,
    config: dict,
    study_context: dict,
    *,
    report_dir: str,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> ReportReviewerResult | None
"""

from __future__ import annotations

import json
import logging
import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm

logger = logging.getLogger(__name__)

# Maximum plain-text characters passed to the LLM (~6000 tokens at ~4 chars/token)
_MAX_CHARS = 24_000


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ReportFlag:
    category: str = ""      # narrative | figures | methods | conclusions
    severity: str = ""      # info | warning | critical
    description: str = ""
    suggestion: str = ""


@dataclass
class ReportReviewerResult(AiResult):
    report_flags: list[ReportFlag] = field(default_factory=list)
    overall_report_quality: str = ""
    review_path: str | None = None


# ---------------------------------------------------------------------------
# HTML stripping
# ---------------------------------------------------------------------------

def _strip_html(text: str) -> str:
    """
    Strip HTML tags from text.

    Uses BeautifulSoup when available; falls back to regex unconditionally
    so the fallback is always exercised and tested regardless of bs4 presence.
    """
    try:
        from bs4 import BeautifulSoup  # type: ignore
        plain = BeautifulSoup(text, "html.parser").get_text(separator=" ")
    except ImportError:
        plain = re.sub(r"<[^>]+>", " ", text)

    # Collapse whitespace produced by either path
    plain = re.sub(r"\s+", " ", plain).strip()
    return plain


# ---------------------------------------------------------------------------
# Markdown writer
# ---------------------------------------------------------------------------

def _write_review_md(result: ReportReviewerResult, report_dir: str) -> str:
    """Write report_review.md and return the path."""
    out_dir = Path(report_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "report_review.md"

    lines: list[str] = ["# Report Review", ""]

    lines += ["## Overall Quality", result.overall_report_quality, ""]

    if result.report_flags:
        lines.append("## Flags")
        lines.append("")
        for flag in result.report_flags:
            category = flag.category.upper() if flag.category else "GENERAL"
            severity = flag.severity.upper() if flag.severity else "INFO"
            lines.append(f"### {category} — {severity}")
            lines.append(f"**Issue:** {flag.description}")
            lines.append(f"**Suggestion:** {flag.suggestion}")
            lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return str(out_path)


# ---------------------------------------------------------------------------
# JSON parser
# ---------------------------------------------------------------------------

def _parse_response(raw: str) -> tuple[list[ReportFlag], str, str, bool]:
    """
    Parse the raw LLM response string into flags, quality paragraph, reasoning.

    Returns (flags, overall_report_quality, reasoning, parse_success).
    On failure returns ([], raw[:500], "", False).
    """
    # Strip optional markdown fences
    cleaned = re.sub(r"^```(?:json)?\s*", "", raw.strip())
    cleaned = re.sub(r"\s*```$", "", cleaned)

    try:
        data = json.loads(cleaned)
    except (json.JSONDecodeError, ValueError):
        logger.warning("report_reviewer: failed to parse LLM response as JSON")
        return [], raw[:500], "", False

    flags: list[ReportFlag] = []
    for item in data.get("report_flags", []):
        if isinstance(item, dict):
            flags.append(ReportFlag(
                category=item.get("category", ""),
                severity=item.get("severity", ""),
                description=item.get("description", ""),
                suggestion=item.get("suggestion", ""),
            ))

    quality = data.get("overall_report_quality", "")
    reasoning = data.get("reasoning", "")
    return flags, quality, reasoning, True


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(
    html_report_path: str,
    config: dict,
    study_context: dict,
    *,
    report_dir: str,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> ReportReviewerResult | None:
    """
    Review the HTML report for scientific quality using the LLM.

    Parameters
    ----------
    html_report_path : str
        Path to the HTML report file to review.
    config : dict
        Pipeline config containing the ``ai`` section.
    study_context : dict
        Parsed study_context.yaml contents.
    report_dir : str
        Output directory for this dataset run (report_review.md written here).
    log_dir : str
        Directory for LLM audit logs.
    runtime_ai : bool
        False → return None immediately (per-step runtime disable).

    Returns
    -------
    ReportReviewerResult or None
    """
    # --- gate ---
    try:
        check_ai_enabled(config, module="report_reviewer", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # --- read HTML ---
    html_path = Path(html_report_path)
    if not html_path.exists():
        logger.warning("report_reviewer: HTML report not found: %s", html_report_path)
        return None

    raw_html = html_path.read_text(encoding="utf-8", errors="replace")

    # --- strip tags and truncate ---
    plain_text = _strip_html(raw_html)
    plain_text = plain_text[:_MAX_CHARS]

    # --- study context fields ---
    tissue = study_context.get("dataset", {}).get("tissue", "unknown")
    disease_context = study_context.get("disease", {}).get("context", None)
    biological_question = study_context.get("biological_question", "")

    # --- LLM call ---
    inputs = {
        "report_text": plain_text,
        "tissue": tissue,
        "disease_context": str(disease_context) if disease_context else "not applicable",
        "biological_question": biological_question or "not specified",
    }

    raw_response = call_llm(
        "report_reviewer",
        inputs,
        config,
        log_dir=log_dir,
        module="report_reviewer",
        runtime_ai=runtime_ai,
    )

    # --- parse ---
    flags, quality, reasoning, parse_success = _parse_response(raw_response)

    # --- read model/provider from config for AiResult fields ---
    ai_cfg = config.get("ai", {})
    model = ai_cfg.get("model", "llama3")
    provider = ai_cfg.get("provider", "ollama")

    result = ReportReviewerResult(
        timestamp=datetime.now(timezone.utc).isoformat(),
        model=model,
        provider=provider,
        skill_name="report_reviewer",
        skill_version="1.0",
        reasoning=reasoning,
        report_flags=flags,
        overall_report_quality=quality,
        review_path=None,
    )

    # --- write markdown ---
    try:
        review_path = _write_review_md(result, report_dir)
        result.review_path = review_path
    except Exception as exc:  # pragma: no cover
        logger.warning("report_reviewer: failed to write review markdown: %s", exc)

    return result
