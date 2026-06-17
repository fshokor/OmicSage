#!/usr/bin/env python3
"""
download_spatial_benchmark.py — OmicSage Phase 7 / Session B

Downloads all spatial transcriptomics benchmark datasets needed for
OmicSage testing and validation.

Usage
-----
    conda activate omicsage
    python scripts/download_spatial_benchmark.py

    # Skip datasets you already have:
    python scripts/download_spatial_benchmark.py --skip kuppe
    python scripts/download_spatial_benchmark.py --skip visium_hd xenium

    # Only validate what is already present (no downloads):
    python scripts/download_spatial_benchmark.py --validate-only

Dataset summary
---------------
Dataset            Size    Source              Requires account?
─────────────────────────────────────────────────────────────────
Kuppe Visium       ~1 GB   Figshare            No  (auto-downloaded)
Kuppe snRNA-seq    ~2 GB   Figshare            No  (auto-downloaded)
Visium HD CRC      ~5 GB   10x Genomics        Yes (manual download)
Xenium Breast      ~8 GB   10x Genomics        Yes (manual download)
─────────────────────────────────────────────────────────────────

The 10x Genomics datasets require a free account registration.
This script will print exact download URLs and instructions,
then validate that the files are in the correct locations.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT   = Path(__file__).resolve().parent.parent
BENCHMARK   = REPO_ROOT / "data" / "benchmark"

# ---------------------------------------------------------------------------
# Dataset registry
# ---------------------------------------------------------------------------

DATASETS = {
    # ------------------------------------------------------------------
    # Kuppe et al. 2022 — auto-downloadable via Figshare
    # ------------------------------------------------------------------
    "kuppe_visium": {
        "label": "Kuppe Visium (human heart, control samples)",
        "dest":  BENCHMARK / "kuppe_visium_human_heart_2022_control.h5ad",
        "url":   "https://figshare.com/ndownloader/files/39347357",
        "type":  "auto",
        "approx_mb": 900,
    },
    "kuppe_snrna": {
        "label": "Kuppe snRNA-seq reference (human heart, control)",
        "dest":  BENCHMARK / "kuppe_snRNA_human_heart_2022_control.h5ad",
        "url":   "https://figshare.com/ndownloader/files/39347573",
        "type":  "auto",
        "approx_mb": 1800,
    },

    # ------------------------------------------------------------------
    # Visium HD — Human Colorectal Cancer (10x Genomics, requires account)
    # ------------------------------------------------------------------
    "visium_hd": {
        "label": "Visium HD Human Colorectal Cancer (10x Genomics)",
        "dest":  BENCHMARK / "visium_hd_colorectal",
        "type":  "manual",
        "approx_mb": 5000,
        "instructions": """\
  1. Go to: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc
  2. Create a free account (or log in) if prompted.
  3. Click "Download" → select "Space Ranger output" → download the .tar.gz file.
  4. Extract it:
       tar -xzf <downloaded_file>.tar.gz
  5. Move the extracted directory to:
       {dest}
     The directory must contain: binned_outputs/ spatial/ and metrics_summary.csv
  6. Re-run this script to validate.""",
        "validate_paths": [
            "binned_outputs",
            "spatial",
        ],
    },

    # ------------------------------------------------------------------
    # Xenium — Human Breast Cancer (10x Genomics, requires account)
    # ------------------------------------------------------------------
    "xenium": {
        "label": "Xenium Human Breast Cancer (10x Genomics)",
        "dest":  BENCHMARK / "xenium_breast",
        "type":  "manual",
        "approx_mb": 8000,
        "instructions": """\
  1. Go to: https://www.10xgenomics.com/datasets/xenium-prime-5k-human-breast-cancer-ffpe
  2. Create a free account (or log in) if prompted.
  3. Click "Download" → select "Xenium Output Bundle" → download the .zip file.
     Minimum files needed (you can skip morphology images to save space):
       • cell_feature_matrix.h5
       • cells.zarr.zip  (or cells.csv.gz for older bundles)
       • experiment.xenium
       • transcripts.parquet  (or transcripts.csv.gz)
  4. Extract the bundle:
       unzip <downloaded_file>.zip -d {dest}
  5. Re-run this script to validate.""",
        "validate_paths": [
            "experiment.xenium",
            "cell_feature_matrix.h5",
            "transcripts.parquet",
        ],
    },
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

RESET  = "\033[0m"
BOLD   = "\033[1m"
GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
CYAN   = "\033[36m"


def _log(msg: str, color: str = RESET):
    print(f"{color}{msg}{RESET}", flush=True)


def _header(msg: str):
    print(f"\n{BOLD}{CYAN}{'─'*60}{RESET}", flush=True)
    print(f"{BOLD}{CYAN}{msg}{RESET}", flush=True)
    print(f"{BOLD}{CYAN}{'─'*60}{RESET}", flush=True)


def _progress_hook(count, block_size, total_size):
    if total_size <= 0:
        print(f"\r  Downloaded {count * block_size / 1_048_576:.1f} MB...", end="", flush=True)
        return
    pct  = min(100, count * block_size * 100 // total_size)
    done = pct // 2
    bar  = "█" * done + "░" * (50 - done)
    mb   = count * block_size / 1_048_576
    tot  = total_size / 1_048_576
    print(f"\r  [{bar}] {pct:3d}%  {mb:.0f}/{tot:.0f} MB", end="", flush=True)


def _download_file(url: str, dest: Path, label: str) -> bool:
    """Download url → dest with a progress bar. Returns True on success."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    try:
        _log(f"  Downloading {label} ...", CYAN)
        _log(f"  URL: {url}")
        urllib.request.urlretrieve(url, tmp, reporthook=_progress_hook)
        print()  # newline after progress bar
        tmp.rename(dest)
        size_mb = dest.stat().st_size / 1_048_576
        _log(f"  ✓  Saved to {dest}  ({size_mb:.0f} MB)", GREEN)
        return True
    except Exception as exc:
        print()
        if tmp.exists():
            tmp.unlink()
        _log(f"  ✗  Download failed: {exc}", RED)
        return False


def _validate_file(dest: Path) -> bool:
    return dest.exists() and dest.stat().st_size > 1_000_000  # > 1 MB


def _validate_dir(dest: Path, required_paths: list[str]) -> tuple[bool, list[str]]:
    if not dest.is_dir():
        return False, required_paths
    missing = [p for p in required_paths if not (dest / p).exists()]
    return len(missing) == 0, missing


# ---------------------------------------------------------------------------
# Per-dataset handlers
# ---------------------------------------------------------------------------

def _handle_auto(key: str, cfg: dict, skip: bool) -> str:
    """Download an auto-downloadable dataset. Returns status string."""
    dest: Path = cfg["dest"]
    label      = cfg["label"]
    url        = cfg["url"]

    if skip:
        _log(f"  Skipped (--skip {key})", YELLOW)
        return "skipped"

    if _validate_file(dest):
        size_mb = dest.stat().st_size / 1_048_576
        _log(f"  ✓  Already present: {dest.name}  ({size_mb:.0f} MB)", GREEN)
        return "ok"

    ok = _download_file(url, dest, label)
    if ok and _validate_file(dest):
        return "ok"
    return "failed"


def _handle_manual(key: str, cfg: dict, skip: bool) -> str:
    """Print manual download instructions. Returns validation status."""
    dest: Path           = cfg["dest"]
    required: list[str]  = cfg.get("validate_paths", [])

    if skip:
        _log(f"  Skipped (--skip {key})", YELLOW)
        return "skipped"

    ok, missing = _validate_dir(dest, required)

    if ok:
        _log(f"  ✓  Already present: {dest}", GREEN)
        for p in required:
            _log(f"    ✓  {p}", GREEN)
        return "ok"

    # Not present — print instructions
    instructions = cfg["instructions"].format(dest=dest)
    _log(f"\n  {BOLD}Manual download required:{RESET}", YELLOW)
    _log(instructions, YELLOW)

    if dest.is_dir() and missing:
        _log(f"\n  Directory exists but missing required files:", RED)
        for p in missing:
            _log(f"    ✗  {p}", RED)

    return "manual_required"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Download OmicSage spatial benchmark datasets.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--skip",
        nargs="+",
        choices=list(DATASETS.keys()),
        default=[],
        metavar="DATASET",
        help="Skip one or more datasets. Choices: %(choices)s",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only check which files are present — no downloads.",
    )
    parser.add_argument(
        "--data-dir",
        default=None,
        help=f"Override benchmark data directory (default: {BENCHMARK})",
    )
    args = parser.parse_args()

    # Allow overriding the benchmark directory
    if args.data_dir:
        bmark = Path(args.data_dir)
        for key, cfg in DATASETS.items():
            cfg["dest"] = bmark / cfg["dest"].name

    skip_set = set(args.skip)

    _header("OmicSage — Spatial Benchmark Dataset Downloader")
    _log(f"Benchmark directory: {BENCHMARK}")
    BENCHMARK.mkdir(parents=True, exist_ok=True)

    results: dict[str, str] = {}

    for key, cfg in DATASETS.items():
        _header(f"{cfg['label']}  (~{cfg['approx_mb']/1000:.0f} GB)")

        if args.validate_only:
            if cfg["type"] == "auto":
                dest = cfg["dest"]
                if _validate_file(dest):
                    size_mb = dest.stat().st_size / 1_048_576
                    _log(f"  ✓  {dest.name}  ({size_mb:.0f} MB)", GREEN)
                    results[key] = "ok"
                else:
                    _log(f"  ✗  Not found: {dest}", RED)
                    results[key] = "missing"
            else:
                ok, missing = _validate_dir(cfg["dest"], cfg.get("validate_paths", []))
                if ok:
                    _log(f"  ✓  {cfg['dest']}", GREEN)
                    results[key] = "ok"
                else:
                    _log(f"  ✗  Not found or incomplete: {cfg['dest']}", RED)
                    if missing:
                        for p in missing:
                            _log(f"    ✗  {p}", RED)
                    results[key] = "missing"
        elif cfg["type"] == "auto":
            results[key] = _handle_auto(key, cfg, key in skip_set)
        else:
            results[key] = _handle_manual(key, cfg, key in skip_set)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    _header("Summary")
    all_ok = True
    for key, status in results.items():
        label = DATASETS[key]["label"]
        if status == "ok":
            _log(f"  ✓  {label}", GREEN)
        elif status == "skipped":
            _log(f"  –  {label}  (skipped)", YELLOW)
        elif status == "manual_required":
            _log(f"  ⚠  {label}  (manual download needed — see instructions above)", YELLOW)
            all_ok = False
        elif status == "failed":
            _log(f"  ✗  {label}  (download failed)", RED)
            all_ok = False
        elif status == "missing":
            _log(f"  ✗  {label}  (not found)", RED)
            all_ok = False

    print()
    if all_ok:
        _log("All datasets ready.", GREEN + BOLD)
    else:
        _log(
            "Some datasets require manual download. "
            "Follow the instructions above, then re-run with --validate-only.",
            YELLOW,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
