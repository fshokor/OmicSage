"""
ai/_llm_client.py
-----------------
Thin BioChatter wrapper for OmicSage.

This is the only file in the AI layer that imports BioChatter directly.
All feature modules call call_llm() — they never touch BioChatter.

Supported providers
-------------------
    claude   →  biochatter.llm_connect.AnthropicConversation
    ollama   →  biochatter.llm_connect.OllamaConversation
    openai   →  biochatter.llm_connect.GptConversation

BioChatter version requirement
-------------------------------
Pinned to biochatter==0.14.2.  Do NOT upgrade without re-running all AI
tests and re-verifying the method names below.  Two confirmed breaking
changes across minor versions:
    - set_system_message()  →  append_system_message()   (renamed in v0.14)
    - base_url is required positional arg for OllamaConversation

Verified method names on v0.14.2
----------------------------------
    append_system_message(text: str)
    query(text: str) → (response: str, token_usage: dict, correction: str)
    set_api_key(api_key: str)          # Claude and OpenAI only

Public API
----------
    call_llm(skill_name, inputs, config, *, log_dir, module,
             runtime_ai) → str

    Returns the raw LLM response string.
    Feature modules are responsible for JSON parsing.

    _build_conversation(config) → BioChatter conversation object
    Exposed for testing (mock injection).

Environment variables
---------------------
    ANTHROPIC_API_KEY   required when provider: claude
    OPENAI_API_KEY      required when provider: openai
    (Ollama needs no key — connects to local server)
"""

from __future__ import annotations

import os
from typing import Any

from ai._audit_log import write_audit_record
from ai._skill_loader import load_skill

# ---------------------------------------------------------------------------
# Valid providers — used in error messages and routing
# ---------------------------------------------------------------------------
_VALID_PROVIDERS = ("claude", "ollama", "openai")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def call_llm(
    skill_name: str,
    inputs: dict[str, Any],
    config: dict[str, Any],
    *,
    log_dir: str = "logs/llm",
    module: str | None = None,
    runtime_ai: bool = True,
) -> str:
    """
    Call the LLM for *skill_name* with *inputs* and return the raw response.

    Parameters
    ----------
    skill_name : str
        Name of the skill YAML in ai/skills/, e.g. "cluster_annotator".
    inputs : dict
        Key-value pairs that fill the skill's user_prompt_template.
    config : dict
        Full pipeline config dict containing an "ai" section.
    log_dir : str, default "logs/llm"
        Directory for audit JSONL files.
    module : str | None
        Module name for the audit log.  Defaults to skill_name if None.
    runtime_ai : bool, default True
        Passed through from the pipeline runner's --ai flag.
        Only used for audit metadata — the gate check is the caller's
        responsibility (done via check_ai_enabled before calling this).

    Returns
    -------
    str
        Raw LLM response string.  May be JSON, plain text, or malformed
        JSON if the model misbehaved.  The caller parses it.

    Raises
    ------
    ValueError
        If config.ai.provider is not one of the valid values.
    KeyError
        If a required environment variable (ANTHROPIC_API_KEY, etc.) is missing.
    Any BioChatter exception
        Network errors, model errors, etc. propagate unchanged so the
        feature module can log them and return None gracefully.
    """
    module_name = module or skill_name
    ai_cfg = config["ai"]
    provider: str = ai_cfg.get("provider", "").lower()
    model: str = ai_cfg.get("model", "")

    # Load skill and fill template
    system_prompt, user_prompt = load_skill(skill_name, **inputs)

    # Build conversation object (injectable for tests via monkeypatching)
    conv = _build_conversation(provider=provider, model=model, config=ai_cfg)

    # Run the query
    conv.append_system_message(system_prompt)
    raw_response, token_usage, _correction = conv.query(user_prompt)

    # Audit — always, no exceptions
    # Derive a compact input_summary: keep scalar values, drop large objects
    input_summary = {
        k: v for k, v in inputs.items()
        if isinstance(v, (str, int, float, bool, type(None)))
    }
    # Count items for list inputs
    for k, v in inputs.items():
        if isinstance(v, list):
            input_summary[f"{k}_count"] = len(v)

    # Extract skill version from the loaded skill file
    skill_version = _get_skill_version(skill_name)

    write_audit_record(
        log_dir=log_dir,
        module=module_name,
        skill_version=skill_version,
        model=model,
        provider=provider,
        input_summary=input_summary,
        token_usage=token_usage if isinstance(token_usage, dict) else None,
        raw_response=raw_response,
        parsed_output=None,       # caller fills this in after parsing
        parse_success=False,      # caller updates audit log after parsing
    )

    return raw_response


def _build_conversation(
    provider: str,
    model: str,
    config: dict[str, Any],
) -> Any:
    """
    Instantiate and return the correct BioChatter conversation object.

    Parameters
    ----------
    provider : str
        "claude" | "ollama" | "openai"
    model : str
        Model identifier string.
    config : dict
        The config["ai"] sub-dict (not the full config).

    Raises
    ------
    ValueError
        Unknown provider.
    KeyError
        Required environment variable missing.
    ImportError
        BioChatter not installed.
    """
    if provider not in _VALID_PROVIDERS:
        raise ValueError(
            f"Unknown AI provider '{provider}'. "
            f"Valid options: {', '.join(_VALID_PROVIDERS)}"
        )

    if provider == "claude":
        from biochatter.llm_connect import AnthropicConversation  # noqa: PLC0415

        api_key = os.environ.get("ANTHROPIC_API_KEY")
        if not api_key:
            raise KeyError(
                "ANTHROPIC_API_KEY environment variable not set. "
                "Export it before running the pipeline with provider: claude."
            )
        conv = AnthropicConversation(
            model_name=model,
            prompts={},
            correct=False,
        )
        conv.set_api_key(api_key=api_key)
        return conv

    if provider == "ollama":
        from biochatter.llm_connect import OllamaConversation  # noqa: PLC0415

        base_url = config.get("ollama_base_url", "http://localhost:11434")
        conv = OllamaConversation(
            base_url=base_url,
            prompts={},
            model_name=model,
            correct=False,
        )
        return conv

    if provider == "openai":
        from biochatter.llm_connect import GptConversation  # noqa: PLC0415

        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            raise KeyError(
                "OPENAI_API_KEY environment variable not set. "
                "Export it before running the pipeline with provider: openai."
            )
        conv = GptConversation(
            model_name=model,
            prompts={},
            correct=False,
        )
        conv.set_api_key(api_key=api_key)
        return conv

    # Unreachable — guard above catches unknown providers
    raise ValueError(f"Unhandled provider: {provider}")  # pragma: no cover


def _get_skill_version(skill_name: str) -> str:
    """
    Read the version field from the skill YAML without re-rendering the template.
    Returns "unknown" if the file or field is missing.
    """
    try:
        import yaml  # noqa: PLC0415
        from pathlib import Path  # noqa: PLC0415

        skill_path = Path("ai/skills") / f"{skill_name}.yaml"
        with skill_path.open("r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
        return str(data.get("version", "unknown"))
    except Exception:  # noqa: BLE001
        return "unknown"
