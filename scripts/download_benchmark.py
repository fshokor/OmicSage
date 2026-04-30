"""
download_benchmark.py
=====================
Downloads and decompresses the GSE194122 NeurIPS 2021 benchmark dataset.
Two files:
  - CITE-seq  : RNA + protein (ADT) from human BMMCs
  - Multiome  : RNA + ATAC from human BMMCs

Usage
-----
    python scripts/download_benchmark.py              # download + unzip both
    python scripts/download_benchmark.py --skip-cite  # multiome only
    python scripts/download_benchmark.py --skip-multiome  # cite only
    python scripts/download_benchmark.py --no-unzip   # keep .gz files

Safe to re-run — skips files that already exist.
"""

import argparse
import gzip
import shutil
import sys
import urllib.request
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent
DEST_DIR   = REPO_ROOT / "data" / "benchmark"

# ── Files to download ──────────────────────────────────────────────────────
BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn"
    "/GSE194122/suppl"
)

FILES = {
    "cite": {
        "gz":   "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz",
        "h5ad": "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad",
        "desc": "CITE-seq  (RNA + protein ADT, ~70k cells)",
    },
    "multiome": {
        "gz":   "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz",
        "h5ad": "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad",
        "desc": "Multiome  (RNA + ATAC, ~69k cells)",
    },
}

# ── Helpers ────────────────────────────────────────────────────────────────

def _progress(block, block_size, total):
    done = block * block_size
    if total > 0:
        pct = min(100, done * 100 // total)
        bar = "█" * (pct // 5) + "░" * (20 - pct // 5)
        mb  = done / 1_048_576
        tot = total / 1_048_576
        print(f"\r  [{bar}] {pct:3d}%  {mb:6.1f} / {tot:.1f} MB",
              end="", flush=True)
    else:
        print(f"\r  {done/1_048_576:6.1f} MB downloaded",
              end="", flush=True)


def download(url: str, dest: Path):
    if dest.exists():
        print(f"  Already exists ({dest.stat().st_size / 1_048_576:.0f} MB) — skipping.")
        return
    print(f"  Downloading:\n    {url}")
    try:
        urllib.request.urlretrieve(url, dest, reporthook=_progress)
        print(f"\n  Saved → {dest.name}  ({dest.stat().st_size / 1_048_576:.0f} MB)")
    except Exception as exc:
        print(f"\n  Download failed: {exc}")
        dest.unlink(missing_ok=True)
        sys.exit(1)


def unzip(gz_path: Path, out_path: Path):
    if out_path.exists():
        print(f"  Already decompressed ({out_path.stat().st_size / 1_048_576:.0f} MB) — skipping.")
        return
    print(f"  Decompressing {gz_path.name} ...")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"  Done → {out_path.name}  ({out_path.stat().st_size / 1_048_576:.0f} MB)")


def verify(h5ad_path: Path):
    """Quick sanity check — file exists and is non-empty."""
    if h5ad_path.exists() and h5ad_path.stat().st_size > 0:
        print(f"  Verified: {h5ad_path.name}")
        return True
    print(f"  MISSING or EMPTY: {h5ad_path.name}")
    return False


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Download GSE194122 NeurIPS 2021 benchmark dataset."
    )
    parser.add_argument("--skip-cite",     action="store_true",
                        help="Skip CITE-seq file")
    parser.add_argument("--skip-multiome", action="store_true",
                        help="Skip Multiome file")
    parser.add_argument("--no-unzip",      action="store_true",
                        help="Keep .gz files, do not decompress")
    parser.add_argument("--keep-gz",       action="store_true",
                        help="Keep .gz after decompression (default: delete)")
    args = parser.parse_args()

    DEST_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 62)
    print("  OmicSage — Benchmark Data Download")
    print("  GSE194122 · NeurIPS 2021 Open Problems in Single-Cell")
    print("=" * 62)
    print(f"  Destination: {DEST_DIR}\n")

    keys = []
    if not args.skip_cite:
        keys.append("cite")
    if not args.skip_multiome:
        keys.append("multiome")

    if not keys:
        print("Nothing to do — both modalities skipped.")
        sys.exit(0)

    all_ok = True

    for key in keys:
        info    = FILES[key]
        gz_path = DEST_DIR / info["gz"]
        h5_path = DEST_DIR / info["h5ad"]

        print(f"\n── {info['desc']} ──")

        # 1. Download .gz
        download(f"{BASE_URL}/{info['gz']}", gz_path)

        # 2. Decompress
        if not args.no_unzip:
            unzip(gz_path, h5_path)

            # 3. Remove .gz unless --keep-gz
            if not args.keep_gz and gz_path.exists():
                gz_path.unlink()
                print(f"  Removed .gz archive to save space.")

            all_ok = all_ok and verify(h5_path)
        else:
            all_ok = all_ok and verify(gz_path)

    print("\n" + "=" * 62)
    if all_ok:
        print("  All files ready.")
        print("\n  Next step:")
        print("    python pipeline/modules/qc/data_report.py \\")
        print("      --input data/benchmark/<file>.h5ad \\")
        print("      --geo GSE194122 \\")
        print("      --output reports/data_intake.html")
    else:
        print("  Some files are missing — re-run this script.")
        sys.exit(1)
    print("=" * 62)


if __name__ == "__main__":
    main()
