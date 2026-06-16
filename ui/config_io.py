"""
OmicSage UI — Config I/O
==========================
Save a config dict to a YAML file and load a YAML config into session state.

Saving: writes to config/runs/<dataset_id>.yaml (matching existing convention).
Loading: parses modality, paths, steps, and params from any OmicSage YAML.
"""
import re
from pathlib import Path
from typing import Any

import yaml


# ── Save ──────────────────────────────────────────────────────────────────────

def save_config(cfg: dict, dataset_id: str) -> str:
    """
    Save cfg to config/runs/<dataset_id>.yaml.
    Creates the directory if it doesn't exist.
    Returns the saved file path.
    """
    out_dir = Path("config/runs")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{dataset_id}.yaml"
    out_path.write_text(
        yaml.dump(cfg, default_flow_style=False, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    return str(out_path)


# ── Load ──────────────────────────────────────────────────────────────────────

def load_config(path: str) -> dict:
    """Load and return a YAML config dict from disk."""
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def parse_config_into_state(cfg: dict, config_path: str) -> dict:
    """
    Parse a loaded YAML config dict into a flat dict suitable for
    injecting into st.session_state.

    Returns a dict with keys matching state.py KEY_* constants' values.
    """
    result: dict[str, Any] = {}

    # ── Detect modality ───────────────────────────────────────────────────────
    # Check dataset.modality, then top-level keys
    ds = cfg.get("dataset", {})
    raw_modality = (
        ds.get("modality", "")
        or ("cite"     if "cite"     in cfg else "")
        or ("multiome" if "multiome" in cfg else "")
        or ("spatial"  if "spatial"  in cfg else "")
        or "scrna"
    )
    modality_map = {
        "scrna":    "scRNA-seq",
        "cite":     "CITE-seq",
        "multiome": "Multiome",
        "spatial":  "Spatial",
        "scRNA-seq": "scRNA-seq",
        "CITE-seq":  "CITE-seq",
        "Multiome":  "Multiome",
        "Spatial":   "Spatial",
    }
    result["modality"]     = modality_map.get(raw_modality, "scRNA-seq")
    result["dataset_name"] = ds.get("name", "") or cfg.get("dataset_id", "")
    result["dataset_id"]   = ds.get("id", "")   or cfg.get("dataset_id", "")
    result["organism"]     = ds.get("organism", "human")

    # ── Paths ─────────────────────────────────────────────────────────────────
    paths = cfg.get("paths", {})
    modality = result["modality"]

    if modality == "scRNA-seq":
        result["data_path"] = paths.get("raw_input", "")
        result["rna_path"]  = ""
    elif modality == "CITE-seq":
        result["data_path"] = paths.get("adt_input", "")
        result["rna_path"]  = paths.get("rna_input", "")
    elif modality == "Multiome":
        result["data_path"] = paths.get("atac_input", "")
        result["rna_path"]  = paths.get("rna_input", "")
    elif modality == "Spatial":
        spatial = cfg.get("spatial", {})
        result["data_path"] = spatial.get("source", paths.get("output_dir", ""))
        result["rna_path"]  = (
            spatial.get("deconvolve", {}).get("ref_path", "")
            or spatial.get("impute", {}).get("sc_reference_path", "")
        )

    result["reports_dir"] = paths.get("reports_dir", "")

    # ── Steps + params ────────────────────────────────────────────────────────
    steps_cfg = cfg.get("steps", {})

    if modality == "Spatial":
        # Spatial uses nested spatial: block, not steps:
        spatial = cfg.get("spatial", {})
        step_keys = ["ingest", "qc", "reduce", "cluster",
                     "deconvolve", "downstream", "impute"]
        selected = []
        step_params = {}
        for s in step_keys:
            block = spatial.get(s, {})
            if isinstance(block, dict):
                # Spatial steps don't have enabled: — presence means enabled
                enabled = block.get("enabled", True)
                if enabled is not False:
                    selected.append(s)
                step_params[s] = {k: v for k, v in block.items() if k != "enabled"}
        result["selected_steps"] = selected
        result["step_params"]    = step_params
    else:
        selected = []
        step_params = {}
        for step, block in steps_cfg.items():
            if isinstance(block, dict) and block.get("enabled", False):
                selected.append(step)
            if isinstance(block, dict):
                step_params[step] = dict(block.get("params", {}))
        result["selected_steps"] = selected
        result["step_params"]    = step_params

    result["config_path"] = config_path
    result["config"]      = cfg

    return result
