"""
adt_reduce.py — ADT dimensionality reduction for OmicSage CITE-seq pipeline.

Takes mdata["adt"].layers["adt_clr"] (CLR-normalized ADT, output of
adt_normalize.py) and performs:
  1. PCA on CLR-normalized ADT values (n_comps capped at min(n_cells-1, n_vars-1, n_comps))
  2. Neighbor graph on ADT PCA embedding (n_pcs=20, per sc-best-practices ch.36)
  3. UMAP (always computed, required)

Outputs written to mdata["adt"]:
  obsm["X_pca_adt"]   — ADT PCA embedding
  obsm["X_umap_adt"]  — ADT UMAP embedding

Naming convention (enforced across OmicSage sessions):
  obsm["X_umap"]      → RNA UMAP  (reduce.py)
  obsm["X_umap_adt"]  → ADT UMAP  (this module)
  obsm["X_umap_wnn"]  → WNN UMAP  (wnn.py, Session 4)

Reference: https://www.sc-best-practices.org/surface_protein/dimensionality_reduction.html
  - sc.pp.pca(svd_solver="arpack")
  - sc.pp.neighbors(n_pcs=20)
  - sc.tl.umap()
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
# Constants (sc-best-practices ch.36 defaults)
# ---------------------------------------------------------------------------
_DEFAULT_N_COMPS = 50       # PCA components computed (variance ratio used to pick n_pcs)
_DEFAULT_N_PCS = 20         # PCs used for neighbor graph (sc-best-practices recommendation)
_DEFAULT_N_NEIGHBORS = 15   # neighbor graph k
_DEFAULT_RANDOM_STATE = 0


def reduce_adt(
    mdata: MuData,
    n_comps: int = _DEFAULT_N_COMPS,
    n_pcs: int = _DEFAULT_N_PCS,
    n_neighbors: int = _DEFAULT_N_NEIGHBORS,
    svd_solver: str = "arpack",
    random_state: int = _DEFAULT_RANDOM_STATE,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Run PCA → neighbor graph → UMAP on CLR-normalized ADT data.

    Parameters
    ----------
    mdata : MuData
        MuData object containing mdata["adt"] with layers["adt_clr"].
    n_comps : int
        Number of PCA components to compute.  Capped at
        min(n_cells - 1, n_vars - 1) automatically.
        Default: 50 (sc-best-practices ch.36).
    n_pcs : int
        Number of PCA components used to build the neighbor graph.
        Must be <= n_comps.  Default: 20 (sc-best-practices ch.36).
    n_neighbors : int
        Number of neighbors for the kNN graph.  Default: 15.
    svd_solver : str
        SVD solver for PCA.  Default: "arpack" (sc-best-practices ch.36).
    random_state : int
        Random seed for reproducibility.  Default: 0.
    inplace : bool
        If True, modify mdata["adt"] in place.
        If False (default), operate on a copy.

    Returns
    -------
    adata : AnnData
        mdata["adt"] (copy or in-place) with:
          obsm["X_pca_adt"]          — ADT PCA embedding  (n_cells × n_comps_actual)
          obsm["X_umap_adt"]         — ADT UMAP embedding (n_cells × 2)
          uns["omicsage_adt_reduce"] — provenance dict
    metrics : dict
        Summary metrics (n_cells, n_vars, n_comps_actual, n_pcs_used,
        variance_explained_total, umap_computed).

    Raises
    ------
    KeyError
        If mdata["adt"] or layers["adt_clr"] is missing.
    ValueError
        If n_pcs > n_comps, or if the ADT matrix is too small for PCA.
    """
    # ------------------------------------------------------------------
    # 1. Input validation
    # ------------------------------------------------------------------
    if "adt" not in mdata.mod:
        raise KeyError("mdata must contain a 'adt' modality. Found: " + str(list(mdata.mod.keys())))

    adt_src = mdata["adt"]

    if "adt_clr" not in adt_src.layers:
        raise KeyError(
            "mdata['adt'].layers['adt_clr'] not found. "
            "Run adt_normalize.py first."
        )

    n_cells, n_vars = adt_src.shape
    if n_cells < 2:
        raise ValueError(f"Need at least 2 cells for PCA, got {n_cells}.")
    if n_vars < 2:
        raise ValueError(f"Need at least 2 ADT features for PCA, got {n_vars}.")

    if n_pcs > n_comps:
        raise ValueError(
            f"n_pcs ({n_pcs}) must be <= n_comps ({n_comps}). "
            "Increase n_comps or decrease n_pcs."
        )

    # ------------------------------------------------------------------
    # 2. Copy or in-place
    # ------------------------------------------------------------------
    adata = adt_src if inplace else adt_src.copy()

    # ------------------------------------------------------------------
    # 3. Set .X to CLR values (scanpy PCA operates on .X)
    # ------------------------------------------------------------------
    adata.X = adata.layers["adt_clr"].copy()

    # ------------------------------------------------------------------
    # 4. Cap n_comps so PCA never requests more components than possible
    # ------------------------------------------------------------------
    max_comps = min(n_cells - 1, n_vars - 1)
    n_comps_actual = min(n_comps, max_comps)

    # Also cap n_pcs if n_comps was reduced below n_pcs
    n_pcs_used = min(n_pcs, n_comps_actual)

    # ------------------------------------------------------------------
    # 5. PCA
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        sc.pp.pca(
            adata,
            n_comps=n_comps_actual,
            svd_solver=svd_solver,
            random_state=random_state,
        )

    # Rename obsm["X_pca"] → obsm["X_pca_adt"] to avoid collision with RNA PCA
    adata.obsm["X_pca_adt"] = adata.obsm.pop("X_pca")

    # Preserve PCA variance info under ADT-specific keys
    if "pca" in adata.uns:
        adata.uns["pca_adt"] = adata.uns.pop("pca")

    # ------------------------------------------------------------------
    # 6. Neighbor graph (uses obsm["X_pca_adt"] via use_rep)
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs_used,
            use_rep="X_pca_adt",
            random_state=random_state,
        )

    # ------------------------------------------------------------------
    # 7. UMAP (always computed — required by OmicSage convention)
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        sc.tl.umap(adata, random_state=random_state)

    # Rename obsm["X_umap"] → obsm["X_umap_adt"]
    adata.obsm["X_umap_adt"] = adata.obsm.pop("X_umap")

    # ------------------------------------------------------------------
    # 8. Variance explained (total across computed components)
    # ------------------------------------------------------------------
    variance_explained_total: Optional[float] = None
    if "pca_adt" in adata.uns and "variance_ratio" in adata.uns["pca_adt"]:
        variance_explained_total = float(
            np.sum(adata.uns["pca_adt"]["variance_ratio"])
        )

    # ------------------------------------------------------------------
    # 9. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "module": "adt_reduce",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "n_comps_requested": n_comps,
            "n_comps_actual": n_comps_actual,
            "n_pcs_used": n_pcs_used,
            "n_neighbors": n_neighbors,
            "svd_solver": svd_solver,
            "random_state": random_state,
        },
        "outputs": {
            "pca_key": "X_pca_adt",
            "umap_key": "X_umap_adt",
        },
    }
    adata.uns["omicsage_adt_reduce"] = provenance

    # ------------------------------------------------------------------
    # 10. Metrics
    # ------------------------------------------------------------------
    metrics: dict = {
        "n_cells": n_cells,
        "n_vars": n_vars,
        "n_comps_actual": n_comps_actual,
        "n_pcs_used": n_pcs_used,
        "n_neighbors": n_neighbors,
        "variance_explained_total": variance_explained_total,
        "umap_computed": True,
        "pca_key": "X_pca_adt",
        "umap_key": "X_umap_adt",
    }

    return adata, metrics
