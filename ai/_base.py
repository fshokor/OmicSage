"""
ai/_base.py
-----------
Base dataclass for all AI feature outputs in OmicSage.

Every feature module (cluster_annotator, pipeline_advisor, etc.) returns
a dataclass that inherits from AiResult.  The base fields carry provenance
and audit metadata so any consumer can trace exactly which model, skill
version, and timestamp produced a given result.

Usage
-----
    from dataclasses import dataclass
    from ai._base import AiResult

    @dataclass
    class ClusterAnnotation(AiResult):
        cell_type: str
        confidence: str   # high | medium | low
        ...
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone


def _utc_now() -> str:
    """Return current UTC time as an ISO-8601 string with Z suffix."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class AiResult:
    """
    Base class for every structured output returned by an OmicSage AI module.

    Fields
    ------
    timestamp : str
        ISO-8601 UTC timestamp of when the LLM call completed.
        Format: "2026-05-12T14:32:01Z"
        Populated automatically by _llm_client.call_llm() — feature modules
        should not set this manually.

    model : str
        Exact model identifier used for this call.
        Examples: "claude-sonnet-4-20250514", "llama3", "gpt-4o"

    provider : str
        Provider name: "claude" | "ollama" | "openai"
        Matches config.ai.provider.

    skill_name : str
        Name of the skill YAML that produced this result.
        Example: "cluster_annotator"

    skill_version : str
        Version string from the skill YAML header.
        Example: "1.0"

    reasoning : str
        The model's free-text reasoning paragraph, carried through from the
        structured JSON response.  Always present; empty string if the skill
        does not request a reasoning field.
    """

    timestamp: str = field(default_factory=_utc_now)
    model: str = ""
    provider: str = ""
    skill_name: str = ""
    skill_version: str = ""
    reasoning: str = ""
