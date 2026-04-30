"""
download_test_data.py
=====================
Downloads GSE166635 (HCC scRNA-seq) into data/test/ for use as
a test fixture for the OmicSage ingestion module.

This dataset provides real 10x MEX files (barcodes/features/matrix)
which are used to test the MTX directory ingestion path in ingest.py.

Usage
-----
    python scripts/download_test_data.py

Data lands in:
    data/test/GSE166635/HCC1/barcodes.tsv.gz
    data/test/GSE166635/HCC1/features.tsv.gz
    data/test/GSE166635/HCC1/matrix.mtx.gz
    data/test/GSE166635/HCC2/barcodes.tsv.gz
    data/test/GSE166635/HCC2/features.tsv.gz
    data/test/GSE166635/HCC2/matrix.mtx.gz

Note: data/test/ is excluded from git via .gitignore (data/ rule).
"""

import sys
import tarfile
import urllib.request
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent
TEST_DIR   = REPO_ROOT / "data" / "test" / "GSE166635"
TAR_FILE   = REPO_ROOT / "data" / "test" / "GSE166635_RAW.tar"

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/"
    "GSE166635/suppl/GSE166635_RAW.tar"
)

EXPECTED = {
    "HCC1": ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"],
    "HCC2": ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"],
}

# ── Helpers ────────────────────────────────────────────────────────────────

def _progress(block, block_size, total):
    done = block * block_size
    if total > 0:
        pct = min(100, done * 100 // total)
        bar = "█" * (pct // 5) + "░" * (20 - pct // 5)
        print(f"\r  [{bar}] {pct:3d}%  {done/1_048_576:.1f}/{total/1_048_576:.1f} MB",
              end="", flush=True)
    else:
        print(f"\r  {done/1_048_576:.1f} MB", end="", flush=True)


def download():
    if TAR_FILE.exists():
        print(f"  Archive already present ({TAR_FILE.stat().st_size/1_048_576:.0f} MB) — skipping.")
        return
    print(f"  Downloading GSE166635_RAW.tar (~204 MB) ...")
    TAR_FILE.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(GEO_URL, TAR_FILE, reporthook=_progress)
        print(f"\n  Saved → {TAR_FILE}")
    except Exception as e:
        print(f"\n  Download failed: {e}")
        TAR_FILE.unlink(missing_ok=True)
        sys.exit(1)


def extract():
    all_present = all(
        (TEST_DIR / sample / fname).exists()
        for sample, fnames in EXPECTED.items()
        for fname in fnames
    )
    if all_present:
        print("  MTX files already extracted — skipping.")
        return

    print(f"  Extracting {TAR_FILE.name} ...")
    with tarfile.open(TAR_FILE, "r") as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            sample = next((s for s in ["HCC1", "HCC2"] if s in name), None)
            if sample is None:
                continue
            if "barcodes" in name.lower():
                dest_name = "barcodes.tsv.gz"
            elif "features" in name.lower() or "genes" in name.lower():
                dest_name = "features.tsv.gz"
            elif "matrix" in name.lower():
                dest_name = "matrix.mtx.gz"
            else:
                continue
            dest = TEST_DIR / sample / dest_name
            if dest.exists():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            print(f"    {name}  →  {sample}/{dest_name}")
            member.name = dest_name
            tar.extract(member, path=TEST_DIR / sample)
    print("  Extraction complete.")


def verify():
    print("\n  Verifying ...")
    all_ok = True
    for sample, fnames in EXPECTED.items():
        for fname in fnames:
            fpath = TEST_DIR / sample / fname
            if fpath.exists() and fpath.stat().st_size > 0:
                print(f"  ✓  {sample}/{fname}  ({fpath.stat().st_size//1024} KB)")
            else:
                print(f"  ✗  {sample}/{fname}  MISSING")
                all_ok = False
    return all_ok


def cleanup():
    if TAR_FILE.exists():
        TAR_FILE.unlink()
        print("  Removed .tar archive.")


# ── Main ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 55)
    print("  OmicSage — Test Data Download (GSE166635)")
    print(f"  Destination: {TEST_DIR}")
    print("=" * 55)

    download()
    extract()
    ok = verify()
    cleanup()

    print("\n" + "=" * 55)
    if ok:
        print("  Test data ready.")
        print(f"\n  MTX directories:")
        print(f"    data/test/GSE166635/HCC1/")
        print(f"    data/test/GSE166635/HCC2/")
        print(f"\n  Run tests:")
        print(f"    pytest tests/test_ingest.py -v")
    else:
        print("  Some files missing — re-run this script.")
        sys.exit(1)
    print("=" * 55)
