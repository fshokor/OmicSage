"""
adt_harmony.py — ADT Harmony batch correction for OmicSage CITE-seq pipeline.

Takes mdata["adt"].obsm["X_pca_adt"] (output of adt_reduce.py) and performs:
  1. Harmony integration on ADT PCA embedding, keyed by a batch column in
     mdata["adt"].obs
  2. Writes corrected embedding to obsm["X_pca_harmony_adt"]
  3. Recomputes neighbor graph using X_pca_harmony_adt
  4. Recomputes UMAP → overwrites obsm["X_umap_adt"] with harmony-corrected UMAP
  5. Returns updated AnnData + metrics dict

Embedding key convention enforced across OmicSage sessions:
  obsm["X_pca_adt"]          — ADT PCA before batch correction  (adt_reduce.py)
  obsm["X_pca_harmony_adt"]  — ADT Harmony-corrected PCA        (this module)
  obsm["X_umap_adt"]         — ADT UMAP (overwritten here with harmony UMAP)
  obsm["X_pca"]              — RNA PCA   (reduce.py  — never touched here)
  obsm["X_pca_harmony"]      — RNA Harmony PCA (future RNA integration)
  obsm["X_umap"]             — RNA UMAP  (reduce.py  — never touched here)

Reference:
  https://www.sc-best-practices.org/surface_protein/batch_correction.html

  Key calls from the reference notebook (ch.38):
    sc.external.pp.harmony_integrate(adata=mdata["prot"], key="donor", random_state=0)
    sc.pp.neighbors(mdata["prot"], n_pcs=20, use_rep="X_pca_harmony", random_state=0)
    sc.tl.umap(mdata["prot"], random_state=0)

  OmicSage naming: use_rep="X_pca_harmony_adt" (not "X_pca_harmony") to avoid
  collision with the RNA harmony embedding (X_pca_harmony) from the RNA pipeline.

API
---
run_harmony_adt(
    mdata,
    batch_key,
    n_pcs=20,
    n_neighbors=15,
    random_state=0,
    inplace=False,
)
→ (AnnData, dict)
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from mudata import MuData


# ---------------------------------------------------------------------------
# Constants (sc-best-practices ch.38 defaults)
# ---------------------------------------------------------------------------
_DEFAULT_N_PCS = 20        # PCs used for post-harmony neighbor graph
_DEFAULT_N_NEIGHBORS = 15  # kNN graph k
_DEFAULT_RANDOM_STATE = 0


def run_harmony_adt(
    mdata: MuData,
    batch_key: str,
    n_pcs: int = _DEFAULT_N_PCS,
    n_neighbors: int = _DEFAULT_N_NEIGHBORS,
    random_state: int = _DEFAULT_RANDOM_STATE,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Harmony batch correction on ADT PCA embedding.

    Reads  ``mdata["adt"].obsm["X_pca_adt"]`` and the batch column
    ``mdata["adt"].obs[batch_key]``, runs Harmony integration, and writes
    the corrected embedding plus a recomputed neighbor graph and UMAP.

    Parameters
    ----------
    mdata : MuData
        MuData object containing ``mdata["adt"]`` with:
          - ``obsm["X_pca_adt"]``  (output of adt_reduce.py)
          - ``obs[batch_key]``     (batch/donor column for Harmony)
    batch_key : str
        Column in ``mdata["adt"].obs`` used as the batch variable for
        Harmony.  Typical value for GSE194122: ``"donor"``.
    n_pcs : int
        Number of Harmony-corrected PCs used to build the post-correction
        neighbor graph.  Default: 20 (sc-best-practices ch.38).
    n_neighbors : int
        Number of neighbors for the post-correction kNN graph.  Default: 15.
    random_state : int
        Random seed for Harmony, neighbors, and UMAP.  Default: 0.
    inplace : bool
        If ``True``, modify ``mdata["adt"]`` in place.
        If ``False`` (default), operate on a copy.

    Returns
    -------
    adata : AnnData
        ``mdata["adt"]`` (copy or in-place) with:
          ``obsm["X_pca_harmony_adt"]``    — Harmony-corrected ADT PCA
          ``obsm["X_umap_adt"]``           — UMAP from harmony embedding
                                             (overwrites uncorrected UMAP)
          ``uns["omicsage_adt_harmony"]``  — provenance dict
    metrics : dict
        Summary metrics — see Notes.

    Raises
    ------
    KeyError
        If ``mdata["adt"]`` is missing, or ``obsm["X_pca_adt"]`` is absent.
    ValueError
        If ``batch_key`` is not found in ``mdata["adt"].obs``, or if
        ``n_pcs`` exceeds the number of PCA components available.

    Notes
    -----
    Metrics dict keys:
      n_cells                 int    number of cells processed
      n_pcs_used              int    PCs used for post-harmony neighbors
      n_neighbors             int    kNN k used
      batch_key               str    batch column name
      n_batches               int    number of unique batch values
      batch_values            list   unique batch values (as strings)
      harmony_key             str    obsm key for corrected embedding
                                     (always "X_pca_harmony_adt")
      umap_key                str    obsm key for harmony UMAP
                                     (always "X_umap_adt")
      umap_computed           bool   always True
      random_state            int    seed used
    """
    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    _validate_inputs(mdata, batch_key, n_pcs)

    # ------------------------------------------------------------------
    # 2. Copy or in-place
    # ------------------------------------------------------------------
    adata = mdata["adt"] if inplace else mdata["adt"].copy()

    # ------------------------------------------------------------------
    # 3. Harmony integration
    #    We call harmonypy directly rather than sc.external.pp.harmony_integrate
    #    to be robust across harmonypy versions.  The scanpy wrapper does
    #    harmony_out.Z_corr.T, which is correct for harmonypy<=0.0.9 (where
    #    Z_corr is PCs×cells) but produces the wrong shape for harmonypy>=2.0.0
    #    (where Z_corr is cells×PCs).  Calling harmonypy directly and storing
    #    the result ourselves avoids this fragility.
    # ------------------------------------------------------------------
    try:
        import harmonypy
    except ImportError as exc:
        raise ImportError(
            "harmonypy is required for ADT batch correction. "
            "Install it with: pip install 'harmonypy==0.0.9'"
        ) from exc

    pca_matrix = adata.obsm["X_pca_adt"].astype(np.float64)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        harmony_out = harmonypy.run_harmony(
            pca_matrix,
            adata.obs,
            batch_key,
            random_state=random_state,
        )

    # Z_corr shape differs by harmonypy version:
    #   <=0.0.9 : (n_pcs, n_cells) — transpose to get (n_cells, n_pcs)
    #   >=2.0.0 : (n_cells, n_pcs) — use directly
    z = harmony_out.Z_corr
    if z.shape[0] != adata.n_obs:
        z = z.T  # harmonypy<=0.0.9 path
    adata.obsm["X_pca_harmony_adt"] = z.astype(np.float32)

    # ------------------------------------------------------------------
    # 4. Recompute neighbor graph on Harmony-corrected embedding
    #    n_pcs capped at actual number of harmony components available
    # ------------------------------------------------------------------
    n_harmony_comps = adata.obsm["X_pca_harmony_adt"].shape[1]
    n_pcs_used = min(n_pcs, n_harmony_comps)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs_used,
            use_rep="X_pca_harmony_adt",
            random_state=random_state,
        )

    # ------------------------------------------------------------------
    # 5. Recompute UMAP → overwrite obsm["X_umap_adt"]
    #    The uncorrected UMAP from adt_reduce.py is intentionally replaced;
    #    all downstream steps (annotation, WNN) use the harmony UMAP.
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.umap(adata, random_state=random_state)

    # Rename obsm["X_umap"] → obsm["X_umap_adt"] (scanpy writes "X_umap")
    adata.obsm["X_umap_adt"] = adata.obsm.pop("X_umap")

    # ------------------------------------------------------------------
    # 6. Metrics
    # ------------------------------------------------------------------
    batch_values = sorted(adata.obs[batch_key].astype(str).unique().tolist())
    metrics: dict = {
        "n_cells": int(adata.n_obs),
        "n_pcs_used": n_pcs_used,
        "n_neighbors": n_neighbors,
        "batch_key": batch_key,
        "n_batches": len(batch_values),
        "batch_values": batch_values,
        "harmony_key": "X_pca_harmony_adt",
        "umap_key": "X_umap_adt",
        "umap_computed": True,
        "random_state": random_state,
    }

    # ------------------------------------------------------------------
    # 7. Provenance
    # ------------------------------------------------------------------
    adata.uns["omicsage_adt_harmony"] = {
        "module": "adt_harmony",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "batch_key": batch_key,
            "n_pcs_requested": n_pcs,
            "n_pcs_used": n_pcs_used,
            "n_neighbors": n_neighbors,
            "random_state": random_state,
        },
        "outputs": {
            "harmony_key": "X_pca_harmony_adt",
            "umap_key": "X_umap_adt",
        },
    }

    return adata, metrics


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_inputs(mdata: MuData, batch_key: str, n_pcs: int) -> None:
    """Raise informative errors for common input mistakes."""
    if "adt" not in mdata.mod:
        raise KeyError(
            "mdata must contain an 'adt' modality. "
            "Found: " + str(list(mdata.mod.keys()))
        )

    adata = mdata["adt"]

    if "X_pca_adt" not in adata.obsm:
        raise KeyError(
            "mdata['adt'].obsm['X_pca_adt'] not found. "
            "Run reduce_adt() before run_harmony_adt()."
        )

    if batch_key not in adata.obs.columns:
        raise ValueError(
            f"batch_key='{batch_key}' not found in mdata['adt'].obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    n_comps_available = adata.obsm["X_pca_adt"].shape[1]
    if n_pcs > n_comps_available:
        raise ValueError(
            f"n_pcs={n_pcs} exceeds the number of PCA components available "
            f"in obsm['X_pca_adt'] ({n_comps_available}). "
            f"Reduce n_pcs or rerun reduce_adt() with more components."
        )
