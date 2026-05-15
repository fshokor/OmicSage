"""
ai/_audit_log.py
----------------
Immutable append-only audit log for every LLM call made by OmicSage.

Design rules
------------
- One JSONL file per module under logs/llm/<module>.jsonl
- One line per call — never overwrite, always append
- The log directory is created automatically if it does not exist
- write_audit_record() never raises — if the write fails, it prints a
  stderr warning and continues.  A logging failure must never crash the
  pipeline or suppress a valid AI result.
- No log rotation — files grow indefinitely.  This is intentional for
  Phase 3: every call is a research artifact.  Rotation can be added later.

Record schema
-------------
Every JSONL line is a JSON object with these keys (all always present):

    timestamp       ISO-8601 UTC  "2026-05-12T14:32:01Z"
    module          str           "cluster_annotator"
    skill_version   str           "1.0"
    model           str           "claude-sonnet-4-20250514"
    provider        str           "claude"
    input_summary   dict          caller-supplied summary of inputs (not full adata)
    prompt_tokens   int | None    from BioChatter token_usage, None if unavailable
    completion_tokens int | None  from BioChatter token_usage, None if unavailable
    raw_response    str           full raw string returned by the LLM
    parsed_output   dict | None   the parsed result, None if parse failed
    parse_success   bool          True if parsed_output is not None

Usage
-----
    from ai._audit_log import write_audit_record

    write_audit_record(
        log_dir="logs/llm",
        module="cluster_annotator",
        skill_version="1.0",
        model="claude-sonnet-4-20250514",
        provider="claude",
        input_summary={"tissue": "liver", "cluster_id": "3", "n_markers": 20},
        token_usage={"prompt_tokens": 412, "completion_tokens": 287},
        raw_response=raw_str,
        parsed_output=result_dict,   # or None if parse failed
        parse_success=True,
    )
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def write_audit_record(
    *,
    log_dir: str | Path,
    module: str,
    skill_version: str,
    model: str,
    provider: str,
    input_summary: dict[str, Any],
    token_usage: dict[str, Any] | None,
    raw_response: str,
    parsed_output: dict[str, Any] | None,
    parse_success: bool,
) -> None:
    """
    Append one audit record to logs/llm/<module>.jsonl.

    Parameters
    ----------
    log_dir : str | Path
        Directory that holds all per-module JSONL files.
        Created automatically if absent.
    module : str
        Name of the calling AI module, e.g. "cluster_annotator".
        Used as the filename stem: <log_dir>/<module>.jsonl
    skill_version : str
        Version from the skill YAML header.
    model : str
        Exact model string used for this call.
    provider : str
        Provider: "claude" | "ollama" | "openai"
    input_summary : dict
        Caller-supplied compact summary of inputs.
        Must NOT contain full AnnData objects — summaries only
        (counts, keys, tissue, etc.).
    token_usage : dict | None
        Token counts from BioChatter's return value.
        Expected keys: "prompt_tokens", "completion_tokens".
        Pass None if BioChatter did not return usage data.
    raw_response : str
        Full raw string returned by the LLM before any parsing.
    parsed_output : dict | None
        The parsed structured result, or None if parsing failed.
    parse_success : bool
        True if parsed_output is a valid dict, False otherwise.

    Returns
    -------
    None
        Always returns — never raises.  Write failures go to stderr.
    """
    try:
        log_path = Path(log_dir) / f"{module}.jsonl"
        log_path.parent.mkdir(parents=True, exist_ok=True)

        # Extract token counts safely
        prompt_tokens: int | None = None
        completion_tokens: int | None = None
        if isinstance(token_usage, dict):
            prompt_tokens = token_usage.get("prompt_tokens")
            completion_tokens = token_usage.get("completion_tokens")

        record: dict[str, Any] = {
            "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "module": module,
            "skill_version": skill_version,
            "model": model,
            "provider": provider,
            "input_summary": input_summary,
            "prompt_tokens": prompt_tokens,
            "completion_tokens": completion_tokens,
            "raw_response": raw_response,
            "parsed_output": parsed_output,
            "parse_success": parse_success,
        }

        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(record, ensure_ascii=False) + "\n")

    except Exception as exc:  # noqa: BLE001
        # Logging must never crash the pipeline.  Print to stderr and continue.
        print(
            f"[omicsage audit_log] WARNING: failed to write audit record "
            f"for module='{module}': {exc}",
            file=sys.stderr,
        )
