"""
ai/_skill_loader.py
===================
Loads a skill YAML file, validates inputs against the skill's declared schema,
fills prompt templates, and returns (system_prompt, user_prompt) ready to send
to the LLM client.

Every AI module calls load_skill() — no raw prompt strings anywhere in Python.

Usage
-----
from ai._skill_loader import load_skill, SkillInputError

system_prompt, user_prompt = load_skill(
    "cluster_annotator",
    tissue="liver",
    disease_context="hepatocellular carcinoma",
    cluster_id="3",
    n_cells=412,
    marker_genes=["CD8A", "GZMB", "PRF1", "PDCD1", "LAG3"],
    prior_knowledge=None,
)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------

class SkillNotFoundError(FileNotFoundError):
    """Raised when the requested skill YAML does not exist."""


class SkillInputError(ValueError):
    """Raised when a required skill input is missing or fails validation."""


class SkillOutputError(ValueError):
    """Raised when the LLM response is missing a required output key."""


# ---------------------------------------------------------------------------
# Skill directory
# ---------------------------------------------------------------------------

_SKILLS_DIR = Path(__file__).parent / "skills"


def _skill_path(skill_name: str) -> Path:
    path = _SKILLS_DIR / f"{skill_name}.yaml"
    if not path.exists():
        available = [p.stem for p in _SKILLS_DIR.glob("*.yaml")]
        raise SkillNotFoundError(
            f"Skill '{skill_name}' not found in {_SKILLS_DIR}. "
            f"Available skills: {available}"
        )
    return path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_skill(skill_name: str, **kwargs: Any) -> tuple[str, str]:
    """
    Load a skill by name, validate inputs, fill templates.

    Parameters
    ----------
    skill_name : str
        Name of the skill (matches filename without .yaml extension).
    **kwargs
        Input values for the skill. Must satisfy the skill's input schema.

    Returns
    -------
    tuple[str, str]
        (system_prompt, user_prompt) ready to send to the LLM.

    Raises
    ------
    SkillNotFoundError
        If the skill YAML does not exist.
    SkillInputError
        If a required input is missing, null, or fails validation.
    """
    skill = _load_yaml(skill_name)
    filled = _validate_and_fill_defaults(skill, kwargs)
    system_prompt = skill["system_prompt"].strip()
    user_prompt = _fill_user_prompt(skill, filled)
    return system_prompt, user_prompt


def load_skill_definition(skill_name: str) -> dict:
    """
    Return the raw skill definition dict.
    Useful for logging the skill version alongside audit records.
    """
    return _load_yaml(skill_name)


def validate_output(skill_name: str, output: dict) -> None:
    """
    Validate a parsed LLM response against the skill's output_schema.

    Parameters
    ----------
    skill_name : str
        Skill to validate against.
    output : dict
        Parsed JSON response from the LLM.

    Raises
    ------
    SkillOutputError
        If a required output key is missing.
    """
    skill = _load_yaml(skill_name)
    schema = skill.get("output_schema", {})
    missing = [key for key in schema if key not in output]
    if missing:
        raise SkillOutputError(
            f"LLM response for skill '{skill_name}' is missing required "
            f"output keys: {missing}. Got keys: {list(output.keys())}"
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_yaml(skill_name: str) -> dict:
    path = _skill_path(skill_name)
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _validate_and_fill_defaults(skill: dict, kwargs: dict) -> dict:
    """
    Validate kwargs against the skill's input schema.
    Fill defaults for optional fields that are not provided.
    Returns a complete dict of filled input values.
    """
    schema = skill.get("inputs", [])
    filled: dict = {}

    for field in schema:
        name = field["name"]
        required = field.get("required", True)
        default = field.get("default", None)
        min_length = field.get("min_length", None)

        if name in kwargs:
            value = kwargs[name]
        elif not required:
            value = default
        else:
            raise SkillInputError(
                f"Skill '{skill['name']}' requires input '{name}' "
                f"({field.get('description', '')})"
            )

        # Validate list min_length
        if min_length is not None and isinstance(value, list):
            if len(value) < min_length:
                raise SkillInputError(
                    f"Skill '{skill['name']}' input '{name}' requires at least "
                    f"{min_length} items, got {len(value)}."
                )

        filled[name] = value

    return filled


def _fill_user_prompt(skill: dict, filled: dict) -> str:
    """
    Fill the user_prompt_template with validated input values.
    Handles list formatting and conditional blocks.
    """
    template: str = skill["user_prompt_template"]

    # Resolve conditional blocks first
    conditional_blocks = skill.get("conditional_blocks", {})
    for block_name, block_def in conditional_blocks.items():
        condition = block_def.get("condition", "")
        default = block_def.get("default", "")
        content_template = block_def.get("content", "")

        # Evaluate simple "X is not null" conditions
        resolved_value = _evaluate_condition(condition, filled)
        if resolved_value:
            # Fill the block content with the same inputs
            block_content = _format_value_in_template(content_template, filled)
        else:
            block_content = default

        filled[block_name] = block_content

    # Format list inputs as numbered or bulleted strings
    format_ready: dict = {}
    for key, value in filled.items():
        format_ready[key] = _format_value(value)

    # Fill template — use a safe approach that doesn't choke on
    # literal braces in the JSON schema example inside the template
    result = _safe_format(template, format_ready)
    return result.strip()


def _evaluate_condition(condition: str, filled: dict) -> bool:
    """
    Evaluate a simple condition string like "prior_knowledge is not null".
    Only supports "X is not null" and "X is null" patterns.
    """
    m = re.match(r"^(\w+)\s+is\s+(not\s+)?null$", condition.strip())
    if not m:
        return False
    field_name = m.group(1)
    negated = m.group(2) is not None  # "not null" → negated=True
    value = filled.get(field_name)
    is_null = value is None or value == "null"
    return negated != is_null  # XOR: "not null" → True when value exists


def _format_value(value: Any) -> str:
    """Format a Python value for insertion into a prompt template."""
    if isinstance(value, list):
        return "\n".join(f"  - {item}" for item in value)
    if value is None:
        return "not specified"
    return str(value)


def _format_value_in_template(template: str, filled: dict) -> str:
    """Fill a sub-template (conditional block content) with filled values."""
    format_ready = {k: _format_value(v) for k, v in filled.items()}
    return _safe_format(template, format_ready)


def _safe_format(template: str, values: dict) -> str:
    """
    Format a template string, leaving literal {{ }} braces (used in the JSON
    schema example in the prompt) untouched.

    Strategy: replace {{ and }} with placeholders, do .format(), restore them.
    """
    placeholder_open  = "__LBRACE__"
    placeholder_close = "__RBRACE__"
    escaped = template.replace("{{", placeholder_open).replace("}}", placeholder_close)

    try:
        filled = escaped.format(**values)
    except KeyError as exc:
        raise SkillInputError(
            f"Prompt template references undefined variable {exc}. "
            f"Available inputs: {list(values.keys())}"
        ) from exc

    result = filled.replace(placeholder_open, "{").replace(placeholder_close, "}")
    return result
