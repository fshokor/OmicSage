"""
ai/_json_utils.py
-----------------
Shared utility for extracting JSON from LLM responses.

Problem solved
--------------
Local models (llama3, mistral, etc.) frequently wrap their JSON output in
prose and markdown fences instead of responding with bare JSON:

    "Here is the result:\n\n```json\n{ ... }\n```\nHope this helps!"

All OmicSage skill parsers call `extract_json_from_response()` before
running json.loads() so they are robust to this pattern.

Extraction strategies (applied in order, first match wins)
-----------------------------------------------------------
1. Bare JSON  — response starts with { or [
2. Fenced block — finds ```(json)? ... ``` anywhere in the response
3. Brace/bracket search — first { to last }, or first [ to last ]

Public API
----------
    extract_json_from_response(raw: str) -> str
        Returns the extracted JSON *string* (not yet parsed).
        Raises ValueError if no JSON block could be found.
"""

from __future__ import annotations

import re


def extract_json_from_response(raw: str) -> str:
    """
    Extract a JSON block from a raw LLM response string.

    Parameters
    ----------
    raw : str
        The raw string returned by the LLM.

    Returns
    -------
    str
        The extracted JSON string, ready for json.loads().

    Raises
    ------
    ValueError
        If the response is empty or contains no extractable JSON.
    """
    text = raw.strip()

    if not text:
        raise ValueError("LLM returned an empty response")

    # ── Strategy 1: response is already bare JSON ──────────────────────────
    if text[0] in ('{', '['):
        return text

    # ── Strategy 2: JSON inside a code fence anywhere in the response ──────
    # Matches ```json or ``` followed by newline, captures up to closing ```
    fence_match = re.search(r'```(?:json)?\s*\n([\s\S]+?)\n\s*```', text)
    if fence_match:
        candidate = fence_match.group(1).strip()
        if candidate:
            return candidate

    # ── Strategy 3: greedy brace extraction — first { to last } ───────────
    first_brace = text.find('{')
    last_brace = text.rfind('}')
    if first_brace >= 0 and last_brace > first_brace:
        return text[first_brace: last_brace + 1]

    # ── Strategy 4: greedy bracket extraction — first [ to last ] ──────────
    first_bracket = text.find('[')
    last_bracket = text.rfind(']')
    if first_bracket >= 0 and last_bracket > first_bracket:
        return text[first_bracket: last_bracket + 1]

    raise ValueError(
        f"No JSON block found in LLM response. First 200 chars: {text[:200]!r}"
    )
