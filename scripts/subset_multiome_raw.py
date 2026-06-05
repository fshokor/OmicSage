"""
subset_multiome_raw.py
======================
Subset the GSE194122 multiome raw AnnData to 6 batches.
Same structure as the input — just fewer cells.

Input  : data/benchmark/GSE194122_multiome_raw_only.h5ad
Output : data/benchmark/GSE194122_multiome_raw_only_6batch.h5ad

Usage
-----
  conda activate omicsage
  python subset_multiome_raw.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_INPUT   = "data/benchmark/GSE194122_multiome_raw_only.h5ad"
DEFAULT_OUTPUT  = "data/benchmark/GSE194122_multiome_raw_only_1batch.h5ad"
DEFAULT_BATCHES = ["s4d8"]#, "s4d1", "s3d10", "s1d2"]


def main(input_path: str, output_path: str, batches: list[str]) -> None:
    import anndata as ad

    print(f"Reading {input_path} ...")
    adata = ad.read_h5ad(input_path)
    print(f"  Full object: {adata.n_obs:,} cells × {adata.n_vars:,} features")

    subset = adata[adata.obs["batch"].isin(batches)].copy()
    print(f"  After subset: {subset.n_obs:,} cells")
    print(subset.obs["batch"].value_counts().to_string())

    print(f"\nWriting {output_path} ...")
    subset.write_h5ad(output_path)
    print("Done.")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",   default=DEFAULT_INPUT)
    parser.add_argument("--output",  default=DEFAULT_OUTPUT)
    parser.add_argument("--batches", nargs="+", default=DEFAULT_BATCHES)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args.input, args.output, args.batches)
