"""
OmicSage — Dataset Preparation Script
scripts/prepare_datasets.py

Strips all processed/computed information from public benchmark datasets,
keeping only raw counts and biological metadata needed for validation.

This must be run ONCE before starting the pipeline. It produces clean
raw-count h5ad files that become the true starting point for OmicSage.

Output files:
  data/benchmark/GSE194122_cite_raw_only.h5ad      (~500 MB)
  data/benchmark/GSE194122_multiome_raw_only.h5ad  (~1.5 GB)

Usage:
    cd ~/OmicSage
    conda activate omicsage
    python scripts/prepare_datasets.py

    # Or process one dataset only:
    python scripts/prepare_datasets.py --dataset cite
    python scripts/prepare_datasets.py --dataset multiome
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BENCHMARK_DIR = Path("data/benchmark")

CITE_INPUT    = BENCHMARK_DIR / "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
MULTIOME_INPUT = BENCHMARK_DIR / "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"

CITE_OUTPUT    = BENCHMARK_DIR / "GSE194122_cite_raw_only.h5ad"
MULTIOME_OUTPUT = BENCHMARK_DIR / "GSE194122_multiome_raw_only.h5ad"


# ---------------------------------------------------------------------------
# Columns to keep per dataset
# ---------------------------------------------------------------------------

# obs columns that are biological metadata — NOT computed by any pipeline step.
# These are kept for:
#   - Validation  : cell_type, GEX_pct_counts_mt
#   - Batch corr  : batch, DonorID, Site, Samplename
#   - Donor info  : everything else

_SHARED_OBS_KEEP = [
    # Ground truth validation
    "cell_type",
    "GEX_pct_counts_mt",        # validate our MT% (r > 0.99)
    # Batch correction
    "batch", "DonorID", "Site", "Samplename",
    # Donor metadata
    "DonorNumber", "DonorAge", "DonorBMI", "DonorBloodType",
    "DonorRace", "Ethnicity", "DonorGender", "DonorSmoker",
    # Other
    "Modality", "is_train",
]

# Multiome has ATAC QC columns — keep raw ATAC metrics for scATAC QC validation
_MULTIOME_EXTRA_OBS_KEEP = [
    "ATAC_nCount_peaks",
    "ATAC_atac_fragments",
    "ATAC_reads_in_peaks_frac",
    "ATAC_blacklist_fraction",
    "ATAC_nucleosome_signal",
]

# What we REMOVE and why:
#   GEX_n_genes_by_counts, GEX_size_factors, GEX_phase  → recomputed by our QC/normalize
#   ADT_* computed columns                               → recomputed by our pipeline
#   All obsm (PCA, UMAP, LSI)                           → recomputed, then compared
#   layers['normalized']                                 → redone by our normalize.py
#   uns except provenance                                → recomputed per step


# ---------------------------------------------------------------------------
# Core stripping function
# ---------------------------------------------------------------------------

def strip_to_raw(
    input_path: Path,
    output_path: Path,
    sample_name: str,
    extra_obs_keep: list[str] | None = None,
) -> ad.AnnData:
    """
    Load a processed public h5ad, strip all computed outputs, save raw version.

    Parameters
    ----------
    input_path:
        Path to the original processed h5ad file.
    output_path:
        Where to write the stripped raw h5ad.
    sample_name:
        Label stored in obs['sample'].
    extra_obs_keep:
        Additional obs columns to keep beyond the shared list.

    Returns
    -------
    Stripped AnnData (also written to output_path).
    """
    logger.info("Loading: %s", input_path.name)
    logger.info("File size: %.1f GB", input_path.stat().st_size / 1e9)

    adata = ad.read_h5ad(input_path)
    logger.info("Loaded: %d cells × %d features", adata.n_obs, adata.n_vars)

    # ── Validate raw counts exist ────────────────────────────────────────
    if "counts" not in adata.layers:
        raise ValueError(
            f"No 'counts' layer found in {input_path.name}. "
            f"Available layers: {list(adata.layers.keys())}"
        )

    raw_counts = adata.layers["counts"]
    logger.info("Raw counts layer found: dtype=%s, sparse=%s",
                raw_counts.dtype, sp.issparse(raw_counts))

    # Verify they are actually integers
    sample = raw_counts.data[:10_000] if sp.issparse(raw_counts) else raw_counts.flatten()[:10_000]
    if not np.all(np.array(sample) == np.floor(np.array(sample))):
        logger.warning("layers['counts'] does not appear to contain integer values — proceed with caution")

    # ── Build obs keep list ──────────────────────────────────────────────
    obs_keep = list(_SHARED_OBS_KEEP)
    if extra_obs_keep:
        obs_keep.extend(extra_obs_keep)

    # Only keep columns that actually exist
    obs_keep_existing = [c for c in obs_keep if c in adata.obs.columns]
    obs_dropped = [c for c in adata.obs.columns if c not in obs_keep_existing and c != "sample"]
    logger.info("Keeping %d obs columns: %s", len(obs_keep_existing), obs_keep_existing)
    logger.info("Dropping %d computed obs columns: %s", len(obs_dropped), obs_dropped)

    # ── Build clean AnnData ──────────────────────────────────────────────
    X = raw_counts.copy()
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.astype(np.float32)

    adata_raw = ad.AnnData(
        X   = X,
        obs = adata.obs[obs_keep_existing].copy(),
        var = adata.var.copy(),        # keep feature_types, gene_id
    )
    adata_raw.obs["sample"] = sample_name

    # Minimal uns — provenance only
    adata_raw.uns["omicsage_source"] = {
        "original_file":  str(input_path.resolve()),
        "sample":         sample_name,
        "stripped":       True,
        "raw_from_layer": "counts",
        "obsm_removed":   list(adata.obsm.keys()),
        "layers_removed": [k for k in adata.layers.keys() if k != "counts"],
    }

    # ── Summary ──────────────────────────────────────────────────────────
    logger.info("─" * 50)
    logger.info("Stripped AnnData: %d cells × %d features", adata_raw.n_obs, adata_raw.n_vars)
    logger.info("X dtype : %s  |  sparse: %s", adata_raw.X.dtype, sp.issparse(adata_raw.X))
    logger.info("obs     : %s", list(adata_raw.obs.columns))
    logger.info("var     : %s", list(adata_raw.var.columns))
    logger.info("obsm    : none (removed: %s)", list(adata.obsm.keys()))
    logger.info("layers  : none (raw counts moved to X)")
    logger.info("─" * 50)

    # ── Save ─────────────────────────────────────────────────────────────
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Saving to: %s", output_path)
    adata_raw.write_h5ad(output_path)
    size_gb = output_path.stat().st_size / 1e9
    logger.info("✓ Saved (%.2f GB)", size_gb)

    return adata_raw


# ---------------------------------------------------------------------------
# Dataset-specific wrappers
# ---------------------------------------------------------------------------

def prepare_cite():
    """Strip GSE194122 CITE-seq processed file."""
    if not CITE_INPUT.exists():
        logger.error("CITE-seq file not found: %s", CITE_INPUT)
        return None

    logger.info("=" * 50)
    logger.info("Preparing CITE-seq dataset")
    logger.info("=" * 50)

    return strip_to_raw(
        input_path   = CITE_INPUT,
        output_path  = CITE_OUTPUT,
        sample_name  = "GSE194122_BMMC_CITE",
        extra_obs_keep = None,   # no extra columns beyond shared list
    )


def prepare_multiome():
    """Strip GSE194122 multiome processed file."""
    if not MULTIOME_INPUT.exists():
        logger.error("Multiome file not found: %s", MULTIOME_INPUT)
        return None

    logger.info("=" * 50)
    logger.info("Preparing Multiome dataset")
    logger.info("=" * 50)

    return strip_to_raw(
        input_path    = MULTIOME_INPUT,
        output_path   = MULTIOME_OUTPUT,
        sample_name   = "GSE194122_BMMC_multiome",
        extra_obs_keep = _MULTIOME_EXTRA_OBS_KEEP,  # keep raw ATAC QC metrics
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Strip processed outputs from public benchmark datasets."
    )
    parser.add_argument(
        "--dataset",
        choices=["cite", "multiome", "all"],
        default="all",
        help="Which dataset to prepare (default: all)",
    )
    args = parser.parse_args()

    if args.dataset in ("cite", "all"):
        prepare_cite()

    if args.dataset in ("multiome", "all"):
        prepare_multiome()

    logger.info("Done. Raw-only files written to %s", BENCHMARK_DIR)


if __name__ == "__main__":
    main()
