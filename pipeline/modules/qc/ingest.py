"""
ingest.py
=========
Data ingestion module for OmicSage.

Auto-detects input format and returns a clean AnnData object
with raw integer counts in adata.X, ready for QC.

Supported formats
-----------------
- .h5ad file     : AnnData format (processed or raw)
- .h5 file       : 10x Genomics HDF5
- directory      : 10x MEX (barcodes.tsv.gz + features.tsv.gz + matrix.mtx.gz)

Key behaviour for processed .h5ad files
----------------------------------------
Many public datasets (e.g. GSE194122) store normalized values in adata.X
and preserve raw counts in a layer (commonly 'counts' or 'spliced').
ingest.py detects this and moves raw counts to adata.X automatically,
storing the normalized matrix in adata.layers['normalized'] for reference.

Usage
-----
    from pipeline.modules.qc.ingest import load_dataset

    adata = load_dataset("data/benchmark/file.h5ad")
    adata = load_dataset("data/benchmark/file.h5")
    adata = load_dataset("data/benchmark/HCC1/")   # MTX directory
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

logger = logging.getLogger(__name__)

# Layers commonly used to store raw counts in processed h5ad files
_RAW_COUNT_LAYER_CANDIDATES = [
    "counts",
    "spliced",
    "raw_counts",
    "UMIs",
    "X_counts",
]


# ── Format detection ────────────────────────────────────────────────────────

def detect_format(path: Path) -> str:
    """
    Detect input format from path.

    Returns one of: 'h5ad', 'h5', 'mtx_dir'
    Raises ValueError for unsupported formats.
    """
    if path.is_dir():
        mtx_files = list(path.glob("matrix.mtx*"))
        if mtx_files:
            return "mtx_dir"
        raise ValueError(
            f"Directory {path} does not look like a 10x MEX folder — "
            f"no matrix.mtx or matrix.mtx.gz found."
        )
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return "h5ad"
    if suffix in (".h5", ".hdf5"):
        return "h5"
    raise ValueError(
        f"Unsupported file format: '{suffix}'. "
        f"Expected .h5ad, .h5, or a 10x MEX directory."
    )


# ── Raw count detection ──────────────────────────────────────────────────────

def _is_integer_matrix(X) -> bool:
    """Return True if matrix contains only integer values."""
    if sp.issparse(X):
        sample = X.data[:10_000] if X.nnz > 10_000 else X.data
    else:
        flat = X.flatten()
        sample = flat[:10_000] if flat.size > 10_000 else flat
    return bool(np.all(sample == np.floor(sample)))


def _find_raw_layer(adata: ad.AnnData) -> Optional[str]:
    """
    Check adata.layers for a layer containing raw integer counts.
    Returns the layer name, or None if not found.
    """
    for candidate in _RAW_COUNT_LAYER_CANDIDATES:
        if candidate in adata.layers:
            if _is_integer_matrix(adata.layers[candidate]):
                logger.info(f"Found raw counts in layer: '{candidate}'")
                return candidate
    # Last resort: check all layers
    for name, layer in adata.layers.items():
        if _is_integer_matrix(layer):
            logger.info(f"Found raw counts in layer: '{name}' (fallback scan)")
            return name
    return None


def _extract_raw_counts(adata: ad.AnnData) -> tuple[ad.AnnData, dict]:
    """
    Ensure adata.X contains raw integer counts.

    Strategy:
    1. If adata.X is already integer → use as-is
    2. If adata.raw exists → pull from adata.raw.X
    3. If a layers contains integers → swap to adata.X
    4. Otherwise → warn and continue (user may know what they're doing)

    Returns (adata, source_info dict)
    """
    source_info = {"raw_layer": None, "x_was_raw": False}

    # Case 1: X is already integer counts
    if _is_integer_matrix(adata.X):
        logger.info("adata.X already contains integer counts — no extraction needed.")
        source_info["x_was_raw"] = True
        return adata, source_info

    # Case 2: adata.raw exists (standard scanpy convention)
    if adata.raw is not None:
        logger.info("Extracting raw counts from adata.raw")
        raw_adata = adata.raw.to_adata()
        # Align vars (raw may have more genes than filtered adata)
        common_vars = adata.var_names.intersection(raw_adata.var_names)
        raw_X = raw_adata[:, common_vars].X
        if _is_integer_matrix(raw_X):
            adata.layers["normalized"] = adata.X.copy()
            adata = adata[:, common_vars].copy()
            adata.X = raw_X
            source_info["raw_layer"] = "adata.raw"
            return adata, source_info

    # Case 3: a layer contains integer counts
    raw_layer = _find_raw_layer(adata)
    if raw_layer is not None:
        logger.info(f"Moving layer '{raw_layer}' → adata.X")
        adata.layers["normalized"] = adata.X.copy()
        adata.X = adata.layers[raw_layer].copy()
        source_info["raw_layer"] = raw_layer
        return adata, source_info

    # Case 4: no raw counts found — warn and continue
    logger.warning(
        "Could not find integer (raw) counts in adata.X, adata.raw, or any layer. "
        "Proceeding with adata.X as-is. QC metrics may be unreliable if data is normalized."
    )
    return adata, source_info


# ── Format-specific loaders ──────────────────────────────────────────────────

def _load_h5ad(path: Path) -> tuple[ad.AnnData, dict]:
    logger.info(f"Loading .h5ad: {path.name}")
    adata = sc.read_h5ad(path)
    adata, source_info = _extract_raw_counts(adata)
    source_info["format"] = "h5ad"
    return adata, source_info


def _load_h5(path: Path) -> tuple[ad.AnnData, dict]:
    logger.info(f"Loading 10x .h5: {path.name}")
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    source_info = {"format": "h5", "raw_layer": None, "x_was_raw": True}
    return adata, source_info


def _load_mtx(path: Path) -> tuple[ad.AnnData, dict]:
    logger.info(f"Loading 10x MEX directory: {path}")
    adata = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)
    source_info = {"format": "mtx_dir", "raw_layer": None, "x_was_raw": True}
    return adata, source_info


# ── Public API ───────────────────────────────────────────────────────────────

def load_dataset(
    path: str | Path,
    sample_name: Optional[str] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Load any supported single-cell dataset and return a clean AnnData
    with raw integer counts in adata.X.

    Parameters
    ----------
    path : str or Path
        Path to .h5ad file, .h5 file, or 10x MEX directory.
    sample_name : str, optional
        Label for this sample, stored in adata.obs['sample'].
        Defaults to the filename stem.
    verbose : bool
        If True, print a summary after loading.

    Returns
    -------
    AnnData
        adata.X             = raw integer counts
        adata.layers['normalized'] = normalized matrix (if input was processed)
        adata.uns['omicsage_source'] = provenance metadata

    Raises
    ------
    FileNotFoundError
        If path does not exist.
    ValueError
        If format is not supported.
    """
    if verbose:
        logging.basicConfig(level=logging.INFO,
                            format="%(levelname)s  %(message)s")

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path not found: {path}")

    fmt = detect_format(path)

    loaders = {
        "h5ad":    _load_h5ad,
        "h5":      _load_h5,
        "mtx_dir": _load_mtx,
    }
    adata, source_info = loaders[fmt](path)

    # Attach sample name
    name = sample_name or path.stem
    adata.obs["sample"] = name

    # Provenance metadata
    adata.uns["omicsage_source"] = {
        "file":       str(path.resolve()),
        "format":     source_info["format"],
        "raw_layer":  source_info.get("raw_layer"),
        "x_was_raw":  source_info.get("x_was_raw", False),
        "sample":     name,
    }

    # Ensure sparse matrix for memory efficiency
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    if verbose:
        _print_summary(adata)

    return adata


def _print_summary(adata: ad.AnnData):
    src = adata.uns["omicsage_source"]
    print("\n" + "─" * 50)
    print(f"  OmicSage · Dataset Loaded")
    print("─" * 50)
    print(f"  Sample   : {src['sample']}")
    print(f"  Format   : {src['format']}")
    print(f"  Cells    : {adata.n_obs:,}")
    print(f"  Features : {adata.n_vars:,}")
    print(f"  Raw from : {'adata.X (already raw)' if src['x_was_raw'] else src['raw_layer']}")
    print(f"  X dtype  : {adata.X.dtype}")
    print(f"  Layers   : {list(adata.layers.keys()) or 'none'}")
    print("─" * 50 + "\n")
