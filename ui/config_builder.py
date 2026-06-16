"""
OmicSage UI — Config Builder
==============================
Converts (modality, data paths, selected steps, step params)
into a valid YAML config dict for the appropriate pipeline runner.

Folder naming logic:
  - dataset.id  = the GEO accession or user-typed ID (e.g. GSE166635)
  - folder name = display name with spaces replaced by underscores,
                  special chars stripped (e.g. "HCC scRNA-seq" → "HCC_scRNA-seq")
  - processed_dir = data/processed/<folder_name>
  - reports_dir   = reports/<folder_name>
"""

import re
import tempfile
from pathlib import Path
from typing import Any

import yaml


def slugify(name: str) -> str:
    """Convert name to lowercase underscore slug — used as fallback ID only."""
    s = re.sub(r"[^\w\s-]", "", name.lower())
    s = re.sub(r"[\s_-]+", "_", s)
    return s.strip("_") or "dataset"


def folder_name(name: str) -> str:
    """
    Convert display name to a filesystem-safe folder name.
    Preserves case and hyphens, only replaces spaces with underscores
    and strips characters that are illegal in folder names.
    e.g. "HCC scRNA-seq (Wang et al. 2025)" → "HCC_scRNA-seq_Wang_et_al_2025"
    """
    s = re.sub(r"[^\w\s\-]", "", name)   # keep letters, digits, _, -, spaces
    s = re.sub(r"\s+", "_", s.strip())   # spaces → underscores
    s = re.sub(r"_+", "_", s)            # collapse multiple underscores
    return s.strip("_") or "dataset"


def write_temp_yaml(cfg: dict) -> str:
    """Write cfg to a named temp file, return path."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", prefix="omicsage_", delete=False
    )
    yaml.dump(cfg, tmp, default_flow_style=False, sort_keys=False, allow_unicode=True)
    tmp.close()
    return tmp.name


# ── scRNA-seq ─────────────────────────────────────────────────────────────────

def build_scrna(dataset_name, data_path, organism,
                selected_steps, step_params, dataset_id="") -> dict:
    did    = dataset_id.strip() if dataset_id and dataset_id.strip() else slugify(dataset_name)
    folder = folder_name(dataset_name)

    cfg: dict[str, Any] = {
        "dataset": {
            "id":       did,
            "name":     dataset_name,
            "modality": "scrna",
            "organism": organism,
        },
        "paths": {
            "raw_input":     data_path,
            "processed_dir": f"data/processed/{folder}",
            "reports_dir":   f"reports/{folder}",
        },
        "steps": {},
    }

    for step in ["qc","normalize","reduce","cluster","annotate",
                 "deg","gsea","harmony","cluster_harmony","pseudobulk"]:
        cfg["steps"][step] = {
            "enabled": step in selected_steps,
            "params":  _clean(step_params.get(step, {})),
        }

    return cfg


# ── CITE-seq ──────────────────────────────────────────────────────────────────

def build_cite(dataset_name, adt_path, rna_path, organism,
               selected_steps, step_params, dataset_id="") -> dict:
    did       = dataset_id.strip() if dataset_id and dataset_id.strip() else slugify(dataset_name)
    folder    = folder_name(dataset_name)
    batch_key = step_params.get("harmony_adt", {}).get("batch_key", "batch")

    cfg: dict[str, Any] = {
        "dataset": {
            "id":       did,
            "name":     dataset_name,
            "modality": "cite",
            "organism": organism,
        },
        "paths": {
            "adt_input":     adt_path,
            "rna_input":     rna_path,
            "processed_dir": f"data/processed/{folder}",
            "reports_dir":   f"reports/{folder}",
        },
        "cite": {"batch_key": batch_key},
        "steps": {},
    }

    for step in ["normalize_adt","doublets","reduce_adt","harmony_adt",
                 "annotate_adt","integration","deg_cite","gsea_cite",
                 "protein_rna_corr","epitope_characterisation"]:
        cfg["steps"][step] = {
            "enabled": step in selected_steps,
            "params":  _clean(step_params.get(step, {})),
        }

    return cfg


# ── Multiome ──────────────────────────────────────────────────────────────────

def build_multiome(dataset_name, atac_path, rna_path, organism,
                   selected_steps, step_params, dataset_id="") -> dict:
    did       = dataset_id.strip() if dataset_id and dataset_id.strip() else slugify(dataset_name)
    folder    = folder_name(dataset_name)
    batch_key = step_params.get("multiome_integration", {}).get("batch_key", "batch")

    cfg: dict[str, Any] = {
        "dataset": {
            "id":   did,
            "name": dataset_name,
        },
        "paths": {
            "atac_input":    atac_path,
            "rna_input":     rna_path,
            "processed_dir": f"data/processed/{folder}",
            "reports_dir":   f"reports/{folder}",
        },
        "multiome": {"batch_key": batch_key},
        "steps": {},
    }

    for step in ["atac_qc","atac_reduce","atac_annotate",
                 "multiome_integration","multiome_deg","multiome_grn"]:
        cfg["steps"][step] = {
            "enabled": step in selected_steps,
            "params":  _clean(step_params.get(step, {})),
        }

    return cfg


# ── Spatial ───────────────────────────────────────────────────────────────────

def build_spatial(dataset_name, data_path, ref_rna_path, organism,
                  selected_steps, step_params, dataset_id="") -> dict:
    did    = dataset_id.strip() if dataset_id and dataset_id.strip() else slugify(dataset_name)
    folder = folder_name(dataset_name)

    ingest_p     = step_params.get("ingest",     {})
    qc_p         = step_params.get("qc",         {})
    reduce_p     = step_params.get("reduce",      {})
    cluster_p    = step_params.get("cluster",     {})
    deconv_p     = step_params.get("deconvolve",  {})
    downstream_p = step_params.get("downstream",  {})
    impute_p     = step_params.get("impute",      {})

    cfg: dict[str, Any] = {
        "dataset_id": did,
        "paths": {
            "output_dir":  f"data/processed/{folder}",
            "reports_dir": f"reports/{folder}",
        },
        "spatial": {
            "source":       data_path,
            "spatial_type": ingest_p.get("spatial_type", "h5ad"),
            "load_images":  ingest_p.get("load_images", True),
            "ingest": {
                "library_key": ingest_p.get("library_key", "sample_name"),
            },
            "qc": _clean({
                "min_counts":   qc_p.get("min_counts",  500),
                "max_counts":   qc_p.get("max_counts",  100000),
                "min_genes":    qc_p.get("min_genes",   200),
                "max_genes":    qc_p.get("max_genes",   10000),
                "max_mt_pct":   qc_p.get("max_mt_pct",  20.0),
                "mt_prefix":    qc_p.get("mt_prefix",   "MT-"),
                "filter_spots": "qc" in selected_steps,
            }),
            "reduce": _clean({
                "n_top_genes":     reduce_p.get("n_top_genes", 3000),
                "n_comps":         reduce_p.get("n_comps",     50),
                "n_neighbors":     reduce_p.get("n_neighbors", 6),
                "target_sum":      reduce_p.get("target_sum",  10000),
                "normalize_total": True,
                "log1p":           True,
                "flavor":          "seurat",
            }) if "reduce" in selected_steps else None,
            "cluster": _clean({
                "resolution":    cluster_p.get("resolution",  0.5),
                "n_neighbors":   cluster_p.get("n_neighbors", 15),
                "n_pcs":         cluster_p.get("n_pcs",       30),
                "random_state":  cluster_p.get("random_state", 0),
                "run_svg":       True,
                "svg_n_genes":   3000,
                "annotation_map": None,
            }) if "cluster" in selected_steps else None,
            "deconvolve": _clean({
                "ref_path":      ref_rna_path,
                "method":        deconv_p.get("method",        "nnls"),
                "cell_type_key": deconv_p.get("cell_type_key", "cell_type_original"),
                "layer_ref":     deconv_p.get("layer_ref",     "counts"),
                "n_jobs":        deconv_p.get("n_jobs",        4),
                "target_sum":    deconv_p.get("target_sum",    10000),
            }) if "deconvolve" in selected_steps else None,
            "downstream": _clean({
                "enabled": "downstream" in selected_steps,
                **{k: downstream_p.get(k, v) for k, v in {
                    "run_region_clustering":   True,
                    "run_celltype_expression": True,
                    "run_co_occurrence":       True,
                    "run_nhood_enrichment":    True,
                    "run_ligrec":              True,
                    "run_svg_gsea":            True,
                    "ligrec_organism":         "human",
                    "n_jobs":                  4,
                }.items()},
            }),
            "impute": _clean({
                "enabled":           "impute" in selected_steps,
                "method":            impute_p.get("method",       "tangram"),
                "n_top_genes":       impute_p.get("n_top_genes",  2000),
                "tangram_mode":      impute_p.get("tangram_mode", "clusters"),
                "device":            impute_p.get("device",       "cpu"),
                "sc_reference_path": ref_rna_path,
            }),
        },
    }

    cfg["spatial"] = {k: v for k, v in cfg["spatial"].items() if v is not None}
    return cfg


# ── Helpers ───────────────────────────────────────────────────────────────────

def _clean(d: dict) -> dict:
    out = {}
    for k, v in d.items():
        if v == "":
            out[k] = None
        elif isinstance(v, dict):
            out[k] = _clean(v)
        else:
            out[k] = v
    return out


def get_reports_dir(cfg: dict, modality: str) -> str:
    return cfg.get("paths", {}).get("reports_dir", "reports")
