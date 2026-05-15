"""
ai/_config_gate.py
------------------
Three-level optionality check for the OmicSage AI layer.

Every AI feature module calls check_ai_enabled() as its very first line.
If AI is disabled at any level, AiDisabledError is raised.  The caller
catches it and returns None — silently, with no log message.

The three levels, checked in order
------------------------------------
Level 1 — global flag:
    config["ai"]["features"] is False  →  disabled for everything

Level 2 — per-module flag:
    config["ai"]["modules"][module_name] is False  →  that module disabled
    If the key is absent, the module defaults to ENABLED.

Level 3 — runtime flag:
    runtime_ai is False  →  disabled regardless of config file
    This corresponds to the --ai flag being absent on run_pipeline.py.
    Passed explicitly by the pipeline runner; defaults to True so that
    notebooks and direct imports work without the flag.

Usage pattern (identical in every feature module)
--------------------------------------------------
    from ai._config_gate import check_ai_enabled, AiDisabledError

    def run(adata, config, runtime_ai: bool = True):
        try:
            check_ai_enabled(config, module="cluster_annotator",
                             runtime_ai=runtime_ai)
        except AiDisabledError:
            return None
        # ... rest of the module — only reached when AI is enabled

Rules
-----
- A disabled module returns None silently. No log. No warning.
- Real errors (network, parse failure, missing key) are NOT caught here.
  They propagate so they can be seen and fixed.
- This module has no imports outside the standard library — it must be
  importable even when biochatter is not installed.
"""

from __future__ import annotations

from typing import Any


class AiDisabledError(Exception):
    """
    Raised by check_ai_enabled() when AI is off at any level.

    Feature modules catch this and return None.  Nothing else should
    catch it — if it propagates past the feature module, that is a bug
    in the feature module, not in the gate.
    """


def check_ai_enabled(
    config: dict[str, Any],
    module: str,
    runtime_ai: bool = True,
) -> None:
    """
    Assert that AI is enabled for *module* given *config* and the runtime flag.

    Parameters
    ----------
    config : dict
        Full pipeline config dict.  Must contain an "ai" top-level key.
        Minimal valid form::

            {"ai": {"features": true, "provider": "claude", "model": "..."}}

    module : str
        Name of the calling module, e.g. "cluster_annotator".
        Must match the key used in config["ai"]["modules"] if that dict
        is present.

    runtime_ai : bool, default True
        Whether the --ai flag was passed on the command line.
        The pipeline runner passes this explicitly.  Notebooks and direct
        imports leave it at the default (True) so they work without the flag.

    Raises
    ------
    AiDisabledError
        If AI is disabled at any of the three levels.

    Returns
    -------
    None
        Returns silently when AI is fully enabled for this module.

    Examples
    --------
    >>> cfg = {"ai": {"features": True, "provider": "claude", "model": "x"}}
    >>> check_ai_enabled(cfg, module="cluster_annotator")   # passes silently

    >>> cfg_off = {"ai": {"features": False}}
    >>> check_ai_enabled(cfg_off, module="cluster_annotator")
    Traceback (most recent call last):
        ...
    ai._config_gate.AiDisabledError: AI disabled globally (ai.features: false)
    """
    # Level 3 — runtime flag (checked first: fastest path, no dict access)
    if not runtime_ai:
        raise AiDisabledError(
            "AI disabled at runtime (--ai flag not passed to run_pipeline.py)"
        )

    # Validate that config has an "ai" section at all
    ai_cfg = config.get("ai")
    if ai_cfg is None:
        raise AiDisabledError(
            "No 'ai' section found in config — AI layer not configured"
        )

    # Level 1 — global features flag
    if not ai_cfg.get("features", True):
        raise AiDisabledError("AI disabled globally (ai.features: false)")

    # Level 2 — per-module flag
    modules_cfg = ai_cfg.get("modules", {})
    if modules_cfg.get(module, True) is False:
        raise AiDisabledError(
            f"AI module '{module}' disabled in config "
            f"(ai.modules.{module}: false)"
        )
