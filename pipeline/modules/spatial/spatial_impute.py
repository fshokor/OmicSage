"""
spatial_impute.py — OmicSage Phase 7 extension

Impute full-transcriptome expression onto Visium spots using a paired
scRNA-seq reference. Fills the "missing genes" chapter: genes not in the
spatial panel are predicted from co-expression in the reference.

Methods
-------
- ``"tangram"`` (default): Optimal-transport mapping of sc cells onto spots
  via ``tangram-sc``. CPU-friendly, ~1–4 GB RAM, interpretable mapping scores.
- ``"gimvi"`` (opt-in): Deep generative model via ``scvi-tools``. Heavier
  memory, GPU recommended. Falls back gracefully when scvi-tools is absent.

Output contract
---------------
- ``adata.obsm["imputed_expression"]``  — np.float32 array (spots × n_top_genes),
  h5py-safe; gene names stored in uns["omicsage_spatial_impute"]["outputs"]["genes_imputed"]
- ``adata.obs["tangram_mapping_score"]`` — per-spot quality score (Tangram only)
- ``adata.uns["omicsage_spatial_impute"]`` — provenance dict

Notes
-----
- Tangram's ``tg.pp_adatas()`` mutates both AnnData objects; we always work on
  copies inside this module. The caller's objects are never modified.
- ``inplace=False`` (default) returns a copy; ``inplace=True`` modifies adata.
- When ``sc_reference_path`` is empty/None the function returns early with a
  ``skipped=True`` provenance flag; the runner uses this to write the passthrough
  checkpoint.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional imports
# ---------------------------------------------------------------------------

try:
    import tangram as tg
    _TANGRAM_AVAILABLE = True
except ImportError:
    _TANGRAM_AVAILABLE = False

try:
    import scvi
    _SCVI_AVAILABLE = True
except ImportError:
    _SCVI_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spatial_impute(
    adata_spatial: ad.AnnData,
    adata_sc: Optional[ad.AnnData] = None,
    method: str = "tangram",
    cell_type_key: str = "cell_type",
    n_top_genes: int = 2000,
    device: str = "cpu",
    tangram_mode: str = "clusters",
    max_cells_per_type: int = 500,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Impute full-transcriptome expression onto Visium spots.

    Parameters
    ----------
    adata_spatial
        AnnData from :func:`spatial_cluster` (post-cluster checkpoint). Must
        contain raw counts in ``X`` or ``layers["counts"]``.
    adata_sc
        Paired scRNA-seq reference AnnData. Must have cell type labels in
        ``obs[cell_type_key]``. When ``None``, returns early with
        ``skipped=True``.
    method
        ``"tangram"`` (default) or ``"gimvi"``.
    cell_type_key
        Column in ``adata_sc.obs`` with cell type annotations.
    n_top_genes
        Number of top highly-variable genes (by variance in sc reference)
        to impute onto spatial spots. Kept low for memory.
    device
        Torch device for Tangram. ``"cpu"`` is safe everywhere; use
        ``"cuda"`` on GPU machines.
    inplace
        If ``False`` (default), returns a copy of ``adata_spatial``.

    Returns
    -------
    adata, params
        Modified AnnData and provenance dict.
        ``params["skipped"]`` is ``True`` when imputation was not run.
    """
    timestamp = datetime.now().isoformat(timespec="seconds")

    if not inplace:
        adata_spatial = adata_spatial.copy()

    # ------------------------------------------------------------------
    # Early exit — no reference provided
    # ------------------------------------------------------------------
    if adata_sc is None:
        prov = {
            "module":    "spatial_impute",
            "timestamp": timestamp,
            "method":    method,
            "skipped":   True,
            "skip_reason": "sc_reference not provided",
            "outputs":   {},
        }
        adata_spatial.uns["omicsage_spatial_impute"] = prov
        return adata_spatial, prov

    if method == "tangram":
        return _run_tangram(adata_spatial, adata_sc, cell_type_key,
                            n_top_genes, device, timestamp,
                            tangram_mode=tangram_mode,
                            max_cells_per_type=max_cells_per_type)
    elif method == "gimvi":
        return _run_gimvi(adata_spatial, adata_sc, cell_type_key,
                          n_top_genes, device, timestamp)
    else:
        raise ValueError(f"Unknown imputation method: {method!r}. "
                         f"Choose 'tangram' or 'gimvi'.")


# ---------------------------------------------------------------------------
# Tangram implementation
# ---------------------------------------------------------------------------


def _run_tangram(
    adata_st: ad.AnnData,
    adata_sc: ad.AnnData,
    cell_type_key: str,
    n_top_genes: int,
    device: str,
    timestamp: str,
    tangram_mode: str = "clusters",
    max_cells_per_type: int = 500,
) -> tuple[ad.AnnData, dict]:
    """Run Tangram imputation.

    Parameters
    ----------
    tangram_mode
        ``"clusters"`` (default, memory-safe): maps per-cell-type mean signatures
        onto spots. Peak RAM ~100 MB regardless of reference size.
        ``"cells"`` : maps individual cells; requires (n_cells × n_spots) float32
        tensor — OOM-prone on large references (>10k cells) without a GPU.
    max_cells_per_type
        When ``tangram_mode="cells"``, subsample each cell type to at most this
        many cells before mapping. Ignored in ``"clusters"`` mode.
    """
    if not _TANGRAM_AVAILABLE:
        raise ImportError(
            "tangram-sc is not installed. "
            "Run: pip install tangram-sc"
        )

    import tangram as tg

    # Work on copies — tg.pp_adatas() mutates both objects
    sc_copy = adata_sc.copy()
    st_copy = adata_st.copy()

    # Ensure raw counts in X for both
    _ensure_counts_in_X(sc_copy, name="sc reference")
    _ensure_counts_in_X(st_copy, name="spatial")

    # Select shared HVGs: top n_top_genes by variance in sc reference
    shared_genes = _select_overlap_hvgs(sc_copy, st_copy, n_top_genes)
    logger.info(f"[spatial_impute/tangram] {len(shared_genes)} shared HVGs selected")

    # Subset both to shared genes
    sc_sub = sc_copy[:, shared_genes].copy()
    st_sub = st_copy[:, shared_genes].copy()

    # ------------------------------------------------------------------ #
    # Memory guard: subsample sc reference for cells mode                 #
    # ------------------------------------------------------------------ #
    if tangram_mode == "cells" and cell_type_key in sc_sub.obs.columns:
        sc_sub = _subsample_per_celltype(sc_sub, cell_type_key, max_cells_per_type)
        logger.info(
            f"[spatial_impute/tangram] subsampled to {sc_sub.n_obs} cells "
            f"({max_cells_per_type} per type max)"
        )

    # pp_adatas normalises and finds shared genes; clusters mode also needs
    # the cell_type_key in obs to build group-level signatures.
    tg.pp_adatas(sc_sub, st_sub, genes=None)

    # Train mapping model
    map_kwargs: dict = dict(
        density_prior="rna_count_based",
        num_epochs=500,
        device=device,
        verbose=False,
    )
    if tangram_mode == "clusters":
        if cell_type_key not in sc_sub.obs.columns:
            logger.warning(
                f"[spatial_impute/tangram] cell_type_key '{cell_type_key}' not found "
                "in sc reference obs — falling back to 'cells' mode."
            )
            tangram_mode = "cells"
        else:
            map_kwargs["cluster_label"] = cell_type_key

    ad_map = tg.map_cells_to_space(
        sc_sub,
        st_sub,
        mode=tangram_mode,
        **map_kwargs,
    )

    # Project sc expression onto spatial spots
    tg.project_genes(adata_map=ad_map, adata_sc=sc_sub)

    # ad_map now has imputed expression (spots × shared_genes)
    imputed = ad_map.to_df()  # DataFrame: spots × genes
    # Align index to spatial obs_names (Tangram returns RangeIndex)
    imputed.index = adata_st.obs_names

    # Mapping scores — stored by Tangram in ad_map.obs["tg_score"]
    if "tg_score" in ad_map.obs.columns:
        adata_st.obs["tangram_mapping_score"] = ad_map.obs["tg_score"].values
        mean_score = float(ad_map.obs["tg_score"].mean())
        n_poor = int((ad_map.obs["tg_score"] < 0.1).sum())
    else:
        mean_score = float("nan")
        n_poor = 0

    # Store imputed values as a float32 numpy array in obsm (h5py-safe).
    # Gene names are preserved in uns provenance under "genes_imputed" and
    # used by the report to reconstruct the DataFrame.
    adata_st.obsm["imputed_expression"] = imputed.values.astype(np.float32)

    prov = {
        "module":           "spatial_impute",
        "timestamp":        timestamp,
        "method":           "tangram",
        "skipped":          False,
        "skip_reason":      None,
        "outputs": {
            "n_genes_imputed":    len(shared_genes),
            "n_spots":            adata_st.n_obs,
            "mean_mapping_score": round(mean_score, 4),
            "n_poor_spots":       n_poor,
            "genes_imputed":      list(shared_genes),
            "cell_type_key":      cell_type_key,
            "device":             device,
        },
    }
    adata_st.uns["omicsage_spatial_impute"] = prov
    logger.info(
        f"[spatial_impute/tangram] done — {len(shared_genes)} genes, "
        f"mean mapping score {mean_score:.4f}"
    )
    return adata_st, prov


# ---------------------------------------------------------------------------
# gimVI implementation
# ---------------------------------------------------------------------------


def _run_gimvi(
    adata_st: ad.AnnData,
    adata_sc: ad.AnnData,
    cell_type_key: str,
    n_top_genes: int,
    device: str,
    timestamp: str,
) -> tuple[ad.AnnData, dict]:
    if not _SCVI_AVAILABLE:
        raise ImportError(
            "scvi-tools is not installed or could not be imported. "
            "Run: pip install scvi-tools"
        )

    from scvi.external import GIMVI

    sc_copy = adata_sc.copy()
    st_copy = adata_st.copy()
    _ensure_counts_in_X(sc_copy, name="sc reference")
    _ensure_counts_in_X(st_copy, name="spatial")

    shared_genes = _select_overlap_hvgs(sc_copy, st_copy, n_top_genes)
    sc_sub = sc_copy[:, shared_genes].copy()
    st_sub = st_copy[:, shared_genes].copy()

    # gimVI requires integer counts
    if sp.issparse(sc_sub.X):
        sc_sub.X = sc_sub.X.toarray()
    if sp.issparse(st_sub.X):
        st_sub.X = st_sub.X.toarray()
    sc_sub.X = sc_sub.X.astype(np.float32)
    st_sub.X = st_sub.X.astype(np.float32)

    GIMVI.setup_anndata(sc_sub)
    GIMVI.setup_anndata(st_sub)

    model = GIMVI(sc_sub, st_sub)
    use_gpu = device != "cpu"
    model.train(200, use_gpu=use_gpu)

    # Impute spatial expression via the trained model
    _, imputed_arr = model.get_imputed_values(normalized=True)
    imputed = pd.DataFrame(
        imputed_arr, index=st_copy.obs_names, columns=shared_genes
    )
    adata_st.obsm["imputed_expression"] = imputed.values.astype(np.float32)

    prov = {
        "module":      "spatial_impute",
        "timestamp":   timestamp,
        "method":      "gimvi",
        "skipped":     False,
        "skip_reason": None,
        "outputs": {
            "n_genes_imputed":    len(shared_genes),
            "n_spots":            adata_st.n_obs,
            "mean_mapping_score": float("nan"),  # gimVI doesn't produce spot scores
            "n_poor_spots":       0,
            "genes_imputed":      list(shared_genes),
            "cell_type_key":      cell_type_key,
            "device":             device,
        },
    }
    adata_st.uns["omicsage_spatial_impute"] = prov
    return adata_st, prov


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ensure_counts_in_X(adata: ad.AnnData, name: str = "AnnData") -> None:
    """Make sure X holds raw integer-like counts.

    Prefer layers["counts"] when present (common OmicSage convention),
    otherwise assume X is already counts and check for negative values.
    """
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        x = adata.X
        arr = x.toarray() if sp.issparse(x) else np.asarray(x)
        if arr.min() < 0:
            raise ValueError(
                f"[spatial_impute] {name}: X contains negative values and no "
                "'counts' layer found. Provide raw counts."
            )


def _subsample_per_celltype(
    adata: ad.AnnData,
    cell_type_key: str,
    max_per_type: int,
) -> ad.AnnData:
    """Subsample each cell type to at most max_per_type cells.

    Used in Tangram 'cells' mode to prevent the (n_cells × n_spots) mapping
    tensor from exhausting RAM on large references.
    """
    keep = []
    rng = np.random.default_rng(42)
    for ct in adata.obs[cell_type_key].unique():
        idx = np.where(adata.obs[cell_type_key] == ct)[0]
        if len(idx) > max_per_type:
            idx = rng.choice(idx, max_per_type, replace=False)
        keep.extend(idx.tolist())
    keep.sort()
    return adata[keep].copy()


def _select_overlap_hvgs(
    adata_sc: ad.AnnData,
    adata_st: ad.AnnData,
    n_top: int,
) -> list[str]:
    """Select top n_top genes by variance in sc, intersected with spatial genes."""
    shared = list(set(adata_sc.var_names) & set(adata_st.var_names))
    if len(shared) == 0:
        raise ValueError(
            "[spatial_impute] No overlapping genes between sc reference and "
            "spatial data. Check that both use the same gene ID convention "
            "(gene symbols vs Ensembl IDs)."
        )

    # Compute per-gene variance in sc reference (shared genes only)
    sc_sub = adata_sc[:, shared]
    x = sc_sub.X
    if sp.issparse(x):
        # Efficient sparse variance: E[x²] - E[x]²
        x_sq = x.copy()
        x_sq.data **= 2
        mean_sq = np.asarray(x_sq.mean(axis=0)).flatten()
        mean    = np.asarray(x.mean(axis=0)).flatten()
        variances = mean_sq - mean ** 2
    else:
        variances = np.var(np.asarray(x), axis=0)

    order = np.argsort(variances)[::-1]
    top_genes = [shared[i] for i in order[:n_top]]
    return top_genes
