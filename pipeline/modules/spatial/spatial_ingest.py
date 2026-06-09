"""
spatial_ingest.py — OmicSage Phase 7, Session 1
Unified spatial transcriptomics data ingestion.

Single entry point for all supported spatial technologies.
Technology is specified via `spatial_type` or auto-detected from the source.

Supported types (current):
  "visium"    — 10x Visium Space Ranger output directory
  "h5ad"      — pre-built AnnData on disk
  "benchmark" — squidpy built-in mouse brain H&E dataset (testing only)

Supported types (future — stubs raise NotImplementedError):
  "visium_hd" — 10x Visium HD (binned_outputs/ directory)
  "xenium"    — 10x Xenium (transcripts.parquet)
  "merfish"   — Vizgen MERSCOPE (cell_by_gene.csv)
  "codex"     — Akoya CODEX / IMC (protein CSV)

Auto-detection fingerprints (when spatial_type="auto"):
  source == "benchmark"                          → benchmark
  source ends with ".h5ad"                       → h5ad
  directory contains spatial/ + .h5 file        → visium
  directory contains binned_outputs/             → visium_hd
  directory contains transcripts.parquet         → xenium
  directory contains cell_by_gene.csv            → merfish
  file ends with .csv (non-spatial dir)          → codex

Output AnnData contract (all implemented types produce):
  obsm["spatial"]                            : coordinates (n_obs, 2)
  uns["spatial"][library_id]["images"]       : tissue images (where available)
  uns["spatial"][library_id]["scalefactors"] : scale factors (where available)
  uns["omicsage_spatial_ingest"]             : provenance dict
    - source, spatial_type, n_obs, n_vars, timestamp, technology_notes
"""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path
from typing import Optional

import anndata as ad
import scanpy as sc

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Supported types registry
# ---------------------------------------------------------------------------

# Maps spatial_type → (loader_function, is_implemented)
# Loaders are defined below; registry is built at module load time.
_LOADER_REGISTRY: dict[str, tuple] = {}   # populated after function defs

# Auto-detection fingerprints: ordered list of (test_fn, spatial_type)
# First match wins.
_AUTO_FINGERPRINTS: list[tuple] = []       # populated after function defs

# Human-readable notes stored in provenance per technology
_TECHNOLOGY_NOTES = {
    "visium":    "10x Visium — spot-based, ~55µm, whole transcriptome, multi-cell resolution",
    "visium_hd": "10x Visium HD — spot-based, ~8µm, whole transcriptome, near single-cell",
    "xenium":    "10x Xenium — imaging-based, single-cell resolution, targeted panel",
    "merfish":   "Vizgen MERSCOPE/MERFISH — imaging-based, single-cell resolution, targeted panel",
    "codex":     "Akoya CODEX / IMC — imaging-based, single-cell resolution, protein markers",
    "h5ad":      "Pre-built AnnData loaded from disk (raw counts preserved; ENSEMBL IDs swapped if var['gene_ids'] present; MT genes stripped)",
    "benchmark": "squidpy built-in mouse brain H&E Visium dataset",
}


# ---------------------------------------------------------------------------
# Image normalisation helper
# ---------------------------------------------------------------------------


def _strip_alpha_from_images(adata: ad.AnnData) -> None:
    """Strip alpha channel from H&E tissue images stored in uns['spatial'].

    Visium h5ad files published on GEO (e.g. Kuppe et al. 2022) store the
    H&E image as an RGBA array with shape (H, W, 4) rather than the RGB
    (H, W, 3) that squidpy's spatial_scatter expects.  When a 4-channel
    array is passed to sq.pl.spatial_scatter it raises a ValueError which
    the report try/except blocks catch silently, producing figures with no
    tissue background.

    This function converts every (H, W, 4) image to (H, W, 3) by dropping
    the alpha channel.  It operates in-place on adata.uns and is idempotent
    (safe to call multiple times; already-RGB images are left unchanged).
    """
    import numpy as np

    for sample_data in adata.uns.get("spatial", {}).values():
        images = sample_data.get("images", {})
        for key, img in list(images.items()):
            if (
                isinstance(img, np.ndarray)
                and img.ndim == 3
                and img.shape[2] == 4
            ):
                images[key] = img[:, :, :3]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spatial_ingest(
    source: str,
    spatial_type: str = "auto",
    counts_file: str = "filtered_feature_bc_matrix.h5",
    library_id: Optional[str] = None,
    library_key: Optional[str] = None,
    load_images: bool = True,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Load spatial transcriptomics data into a standard AnnData.

    Single entry point for all supported spatial technologies.
    The technology is selected via *spatial_type* or auto-detected from
    the structure of *source*.

    Parameters
    ----------
    source
        One of:
        - ``"benchmark"`` — squidpy built-in mouse brain H&E dataset
        - path to a ``.h5ad`` file
        - path to a Space Ranger / Xenium / Vizgen output directory
    spatial_type
        Technology type. One of:
        ``"auto"`` (default), ``"visium"``, ``"visium_hd"``,
        ``"xenium"``, ``"merfish"``, ``"codex"``, ``"h5ad"``,
        ``"benchmark"``.
        With ``"auto"`` the type is inferred from *source*.
    counts_file
        Count matrix filename (Visium only).
        Default: ``"filtered_feature_bc_matrix.h5"``.
    library_id
        Identifier stored in ``uns["spatial"]`` (Visium only).
        When ``None``, inferred from the directory name.
    library_key
        Column in ``adata.obs`` that maps each spot to its library/sample ID.
        Required by ``sq.pl.spatial_scatter`` when multiple samples are merged
        into one AnnData (i.e. ``uns["spatial"]`` has more than one key).
        When ``None`` (default), the value is auto-detected at ingest time by
        matching ``uns["spatial"]`` keys against obs columns, and stored in
        provenance so report generators can use it without re-detecting.
        Set explicitly when auto-detection fails (e.g. sample names were
        renamed during concatenation).
    load_images
        Whether to load tissue images (Visium only).
    inplace
        Ignored — always returns a new object. Present for API consistency.

    Returns
    -------
    adata
        AnnData with ``obsm["spatial"]`` and
        ``uns["omicsage_spatial_ingest"]``.
    params
        Provenance dictionary.

    Raises
    ------
    ImportError
        If squidpy is not installed.
    ValueError
        If the source type cannot be determined or is unsupported.
    NotImplementedError
        If *spatial_type* is recognised but not yet implemented
        (e.g. ``"xenium"``, ``"merfish"``, ``"codex"``).
    """
    _check_squidpy()

    resolved_type = _resolve_spatial_type(source, spatial_type)
    loader, is_implemented = _LOADER_REGISTRY.get(
        resolved_type, (None, False)
    )

    if loader is None:
        raise ValueError(
            f"Unknown spatial_type={resolved_type!r}. "
            f"Supported: {list(_LOADER_REGISTRY)}"
        )
    if not is_implemented:
        raise NotImplementedError(
            f"spatial_type={resolved_type!r} is not yet implemented. "
            f"It is planned for a future OmicSage phase. "
            f"To load this data manually, read it into an AnnData with "
            f"obsm['spatial'] set, save as .h5ad, and use spatial_type='h5ad'."
        )

    adata, effective_library_id, source_repr = loader(
        source,
        counts_file=counts_file,
        library_id=library_id,
        load_images=load_images,
    )

    _validate_spatial_adata(adata, source_repr)

    # Resolve library_key: use explicit value if given, otherwise auto-detect.
    # Stored in provenance so report generators read it directly instead of
    # re-running detection on every report call.
    resolved_library_key = library_key or _detect_library_key(adata)

    params = {
        "source": source_repr,
        "spatial_type": resolved_type,
        "technology_notes": _TECHNOLOGY_NOTES.get(resolved_type, ""),
        "counts_file": counts_file,
        "library_id": effective_library_id,
        "library_key": resolved_library_key,
        "load_images": load_images,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "timestamp": datetime.now().isoformat(),
    }

    adata.uns["omicsage_spatial_ingest"] = params
    return adata, params


def list_supported_types() -> dict[str, str]:
    """Return all spatial types and their implementation status.

    Returns
    -------
    dict mapping spatial_type → "implemented" | "planned"
    """
    return {
        k: ("implemented" if v[1] else "planned")
        for k, v in _LOADER_REGISTRY.items()
    }


def _detect_library_key(adata: ad.AnnData) -> Optional[str]:
    """Auto-detect the obs column that maps spots to their library/sample ID.

    Matches the keys of ``uns["spatial"]`` against obs column values.
    Returns ``None`` for single-sample data (no ``library_key`` needed by
    squidpy) or when no matching column is found.

    The result is stored in ``uns["omicsage_spatial_ingest"]["library_key"]``
    at ingest time so downstream report generators can read it directly
    rather than re-running detection.

    Priority order: common column names checked first for speed, then a full
    scan of all string/category columns as fallback.
    """
    library_ids = list(adata.uns.get("spatial", {}).keys())
    if len(library_ids) <= 1:
        return None  # single-sample — squidpy does not need library_key
    id_set = set(library_ids)
    for candidate in ("library_id", "sample", "patient", "donor_id", "batch", "slide"):
        if candidate in adata.obs.columns:
            if id_set.issubset(set(adata.obs[candidate].astype(str).unique())):
                return candidate
    # Full scan fallback
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name in ("object", "category"):
            if id_set.issubset(set(adata.obs[col].astype(str).unique())):
                return col
    return None


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------


def _resolve_spatial_type(source: str, spatial_type: str) -> str:
    """Return the effective spatial type, either explicit or auto-detected."""
    if spatial_type != "auto":
        return spatial_type

    src = str(source)

    for test_fn, detected_type in _AUTO_FINGERPRINTS:
        if test_fn(src):
            return detected_type

    raise ValueError(
        f"Cannot auto-detect spatial_type for source={source!r}.\n"
        f"Recognised fingerprints:\n"
        f"  'benchmark'            — pass source='benchmark'\n"
        f"  Visium                 — directory with spatial/ subfolder\n"
        f"  Visium HD              — directory with binned_outputs/ subfolder\n"
        f"  Xenium                 — directory with transcripts.parquet\n"
        f"  MERFISH/Vizgen         — directory with cell_by_gene.csv\n"
        f"  CODEX/IMC              — .csv file with protein markers\n"
        f"  Pre-built AnnData      — path ending in .h5ad\n"
        f"Pass spatial_type explicitly to override auto-detection."
    )


def _is_benchmark(src: str) -> bool:
    return src == "benchmark"

def _is_h5ad(src: str) -> bool:
    return src.endswith(".h5ad")

def _is_visium(src: str) -> bool:
    return (
        os.path.isdir(src)
        and os.path.isdir(os.path.join(src, "spatial"))
    )

def _is_visium_hd(src: str) -> bool:
    return (
        os.path.isdir(src)
        and os.path.isdir(os.path.join(src, "binned_outputs"))
    )

def _is_xenium(src: str) -> bool:
    return (
        os.path.isdir(src)
        and os.path.isfile(os.path.join(src, "transcripts.parquet"))
    )

def _is_merfish(src: str) -> bool:
    return (
        os.path.isdir(src)
        and os.path.isfile(os.path.join(src, "cell_by_gene.csv"))
    )

def _is_codex(src: str) -> bool:
    return (
        os.path.isfile(src)
        and src.endswith(".csv")
        and not os.path.isdir(src)
    )


# ---------------------------------------------------------------------------
# Loader functions
# Signature: (source, *, counts_file, library_id, load_images)
#            → (AnnData, effective_library_id, source_repr)
# ---------------------------------------------------------------------------


def _load_benchmark(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    adata = sq.datasets.visium_hne_adata()
    return adata, library_id or "benchmark", "squidpy:visium_hne_adata"


def _load_h5ad(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    path = str(source)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"h5ad file not found: {path!r}")
    adata = sc.read_h5ad(path)

    # ------------------------------------------------------------------ #
    # Normalise to the standard spatial AnnData contract
    # ------------------------------------------------------------------ #

    # 1. Preserve raw counts in layers["counts"].
    #    For Kuppe-style Visium h5ad files, X contains raw counts and
    #    layers is empty — save X before any downstream normalization.
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # 2. Swap var_names to ENSEMBL gene IDs when var["gene_ids"] is present.
    #    cell2location requires ENSEMBL IDs to map between spatial and
    #    reference data.  Space Ranger h5ad files store gene symbols as
    #    var_names and ENSEMBL IDs in var["gene_ids"].
    if "gene_ids" in adata.var.columns:
        adata.var["feature_name"] = adata.var_names.copy()
        adata.var_names = adata.var["gene_ids"].astype(str)
        adata.var_names_make_unique()

    # 3. Strip mitochondrial genes into obsm["MT"].
    #    cell2location recommends removing MT genes before deconvolution
    #    as they represent technical artefacts rather than cell abundance.
    #    MT genes are identified by the "MT-" prefix on the feature_name
    #    column (populated in step 2) or var_names as a fallback.
    if "feature_name" in adata.var.columns:
        # pandas Series.str.startswith → Series → .values → numpy array
        mt_mask = adata.var["feature_name"].str.startswith("MT-").values
    else:
        # pandas Index.str.startswith returns numpy array directly — no .values needed
        mt_mask = adata.var_names.str.startswith("MT-")

    if mt_mask.sum() > 0:
        import scipy.sparse as sp
        mt_matrix = adata[:, mt_mask].X
        if sp.issparse(mt_matrix):
            mt_matrix = mt_matrix.toarray()
        adata.obsm["MT"] = mt_matrix
        adata = adata[:, ~mt_mask].copy()

    # 4. Strip alpha channel from H&E images (RGBA → RGB).
    #    Some published h5ad files store images as (H, W, 4); squidpy
    #    spatial_scatter requires (H, W, 3).
    _strip_alpha_from_images(adata)

    return adata, library_id or "custom", path


def _load_visium(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    path = str(source)
    if not os.path.isdir(path):
        raise NotADirectoryError(
            f"Visium Space Ranger directory not found: {path!r}"
        )
    kwargs: dict = dict(counts_file=counts_file, load_images=load_images)
    if library_id is not None:
        kwargs["library_id"] = library_id

    adata = sq.read.visium(path, **kwargs)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    _strip_alpha_from_images(adata)
    return adata, library_id or Path(path).name, path


def _load_visium_hd(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    raise NotImplementedError(
        "Visium HD support is planned for a future OmicSage phase. "
        "To load manually: use squidpy.read.visium_hd() (squidpy >= 1.5), "
        "then save as .h5ad and reload with spatial_type='h5ad'."
    )


def _load_xenium(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    raise NotImplementedError(
        "Xenium support is planned for a future OmicSage phase. "
        "To load manually: use squidpy.read.xenium(), "
        "then save as .h5ad and reload with spatial_type='h5ad'."
    )


def _load_merfish(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    raise NotImplementedError(
        "MERFISH/Vizgen support is planned for a future OmicSage phase. "
        "To load manually: use squidpy.read.vizgen(), "
        "then save as .h5ad and reload with spatial_type='h5ad'."
    )


def _load_codex(
    source: str, *, counts_file, library_id, load_images
) -> tuple[ad.AnnData, str, str]:
    raise NotImplementedError(
        "CODEX/IMC support is planned for a future OmicSage phase. "
        "To load manually: read the CSV into an AnnData with obsm['spatial'] "
        "set, save as .h5ad and reload with spatial_type='h5ad'."
    )


# ---------------------------------------------------------------------------
# Registry + fingerprint initialisation (after all functions are defined)
# ---------------------------------------------------------------------------

_LOADER_REGISTRY = {
    "benchmark": (_load_benchmark, True),
    "h5ad":      (_load_h5ad,      True),
    "visium":    (_load_visium,    True),
    "visium_hd": (_load_visium_hd, False),   # planned
    "xenium":    (_load_xenium,    False),   # planned
    "merfish":   (_load_merfish,   False),   # planned
    "codex":     (_load_codex,     False),   # planned
}

# Order matters — more specific fingerprints before generic ones
_AUTO_FINGERPRINTS = [
    (_is_benchmark, "benchmark"),
    (_is_h5ad,      "h5ad"),
    (_is_xenium,    "xenium"),    # check before visium (both are dirs)
    (_is_merfish,   "merfish"),   # check before visium
    (_is_visium_hd, "visium_hd"), # check before visium
    (_is_visium,    "visium"),
    (_is_codex,     "codex"),
]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _check_squidpy() -> None:
    if not _SQUIDPY_AVAILABLE:
        raise ImportError(
            "squidpy is required for spatial analysis. "
            "Install with: pip install squidpy"
        )


def _validate_spatial_adata(adata: ad.AnnData, source_repr: str) -> None:
    """Verify the loaded AnnData satisfies the spatial contract."""
    if "spatial" not in adata.obsm:
        raise ValueError(
            f"Loaded AnnData from {source_repr!r} is missing obsm['spatial']. "
            "Spatial coordinates are required for all downstream spatial modules."
        )
    if "spatial" not in adata.uns:
        raise ValueError(
            f"Loaded AnnData from {source_repr!r} is missing uns['spatial']. "
            "Scale factors / image metadata are required for spatial QC plots."
        )
    coords = adata.obsm["spatial"]
    if coords.shape[1] != 2:
        raise ValueError(
            f"obsm['spatial'] must have shape (n_obs, 2), got {coords.shape}."
        )
