"""
spatial_deconvolve.py — OmicSage Phase 7, Session 5 (multi-method)

Cell type deconvolution of Visium spots. Supports two methods plus a
graceful skip mode.

Methods
-------
- ``"nnls"`` (default): Non-negative least squares. Pure Python, ~200 MB
  peak RAM regardless of dataset size. No extra dependencies beyond
  scipy/sklearn/joblib. Returns *proportions* (rows sum to 1).
- ``"cell2location"``: Bayesian deconvolution. Heavy memory (3-20 GB)
  but gold-standard quality. Returns absolute cell-type abundances.
- ``"none"``: Skip deconvolution gracefully (provenance flag only).

Memory control for cell2location
--------------------------------
``per_sample=True`` runs the spatial fit one library at a time, looping
over ``library_key``. The reference is fit once (it's the same across
libraries). For ~4-sample datasets this drops peak RAM by roughly Nx
(e.g. 16 GB → 4 GB). For NNLS this flag is a no-op (already memory-light).

Output contract (identical for both methods so downstream code is method-agnostic)
---------------------------------------------------------------------------------
- ``obsm["cell_type_abundances"]``    — (n_spots, n_cell_types) canonical name
- ``obsm["q05_cell_abundance_w_sf"]`` — alias of the above (back-compat with code
                                        that hard-codes the cell2location name)
- ``obs[<cell_type>]``                — per-cell-type abundance/proportion column
- ``obs["dominant_cell_type"]``       — categorical argmax across cell types
- ``uns["omicsage_spatial_deconvolve"]`` — provenance dict
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

try:
    import cell2location as c2l
    _C2L_AVAILABLE = True
except ImportError:
    _C2L_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spatial_deconvolve(
    adata: ad.AnnData,
    ref_adata: Optional[ad.AnnData] = None,
    method: str = "nnls",
    per_sample: bool = False,
    library_key: Optional[str] = None,
    cell_type_key: str = "cell_type_original",
    layer_ref: str = "counts",
    # NNLS-specific
    n_jobs: int = 4,
    target_sum: float = 1e4,
    # cell2location-specific
    batch_key_ref: Optional[str] = "donor_id",
    batch_key_st: Optional[str] = "patient",
    covariate_keys: Optional[list] = None,
    N_cells_per_location: int = 8,
    detection_alpha: int = 20,
    max_epochs_ref: int = 250,
    max_epochs_st: int = 30000,
    batch_size_ref: int = 2500,
    batch_size_st: Optional[int] = None,
    num_samples_posterior: int = 1000,
    # cell2location gene-selection thresholds
    cell_count_cutoff: int = 5,
    cell_percentage_cutoff2: float = 0.03,
    nonz_mean_cutoff: float = 1.12,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Deconvolve Visium spots into cell type abundances.

    Parameters
    ----------
    adata
        AnnData produced by :func:`spatial_cluster`. Must contain
        ``layers["counts"]`` with raw integer counts.
    ref_adata
        Paired scRNA-seq reference. When ``None``, deconvolution is skipped.
    method
        ``"nnls"`` (default), ``"cell2location"``, or ``"none"``.
    per_sample
        Run cell2location one library at a time. Cuts peak RAM by ~Nx for
        N-sample data. No effect on NNLS.
    library_key
        Column in ``adata.obs`` identifying samples. Required when
        ``per_sample=True``. Auto-resolved from
        ``uns["omicsage_spatial_ingest"]["library_key"]`` when ``None``.
    cell_type_key
        Column in ``ref_adata.obs`` with cell type labels.
    layer_ref
        Layer in ``ref_adata`` with raw counts.
    n_jobs
        Parallelism for NNLS solve. Default 4.
    target_sum
        Library-size normalization target for NNLS. Default 1e4.
    batch_key_ref, batch_key_st, covariate_keys
        cell2location batch correction columns. ``None`` disables.
    N_cells_per_location
        cell2location prior on cells per spot. Default 8 (heart tissue).
    detection_alpha, max_epochs_ref, max_epochs_st, batch_size_ref
        cell2location hyperparameters.
    batch_size_st
        Spatial training batch size. Set to e.g. 2048 to reduce memory.
        Default ``None`` = use full data.
    num_samples_posterior
        Posterior samples drawn at the end of cell2location. Reducing from
        1000 to 200 roughly halves ``export_posterior`` memory.
    cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff
        cell2location gene-selection thresholds.
    inplace
        Modify *adata* in place.

    Returns
    -------
    adata
        AnnData with deconvolution outputs.
    params
        Provenance dictionary.
    """
    _validate_input(adata, ref_adata, cell_type_key, layer_ref, method)
    adata = adata if inplace else adata.copy()

    # ------------------------------------------------------------------ #
    # Graceful skip
    # ------------------------------------------------------------------ #
    if ref_adata is None or method == "none":
        reason = (
            "method='none' — deconvolution explicitly disabled."
            if method == "none"
            else "ref_adata=None — no scRNA-seq reference provided."
        )
        params = _skip_provenance(method, reason)
        adata.uns["omicsage_spatial_deconvolve"] = params
        return adata, params

    # ------------------------------------------------------------------ #
    # Resolve library_key for per-sample mode
    # ------------------------------------------------------------------ #
    if per_sample:
        if library_key is None:
            library_key = (
                adata.uns.get("omicsage_spatial_ingest", {}).get("library_key")
            )
        if library_key is None or library_key not in adata.obs.columns:
            logger.warning(
                "per_sample=True but no usable library_key found. "
                "Falling back to single-pass deconvolution."
            )
            per_sample = False

    # ------------------------------------------------------------------ #
    # Dispatch
    # ------------------------------------------------------------------ #
    common = {
        "method": method,
        "per_sample": per_sample,
        "library_key": library_key,
        "cell_type_key": cell_type_key,
        "layer_ref": layer_ref,
    }

    if method == "nnls":
        proportions, cell_type_names, n_shared = _deconvolve_nnls(
            adata, ref_adata,
            cell_type_key=cell_type_key,
            layer_ref=layer_ref,
            n_jobs=n_jobs,
            target_sum=target_sum,
        )
        method_params = {"n_jobs": n_jobs, "target_sum": target_sum}

    elif method == "cell2location":
        _check_c2l()
        c2l_kwargs = dict(
            cell_type_key=cell_type_key,
            layer_ref=layer_ref,
            batch_key_ref=batch_key_ref,
            covariate_keys=covariate_keys,
            N_cells_per_location=N_cells_per_location,
            detection_alpha=detection_alpha,
            max_epochs_ref=max_epochs_ref,
            max_epochs_st=max_epochs_st,
            batch_size_ref=batch_size_ref,
            batch_size_st=batch_size_st,
            num_samples_posterior=num_samples_posterior,
            cell_count_cutoff=cell_count_cutoff,
            cell_percentage_cutoff2=cell_percentage_cutoff2,
            nonz_mean_cutoff=nonz_mean_cutoff,
        )
        if per_sample:
            proportions, cell_type_names, n_shared = _deconvolve_c2l_per_sample(
                adata, ref_adata, library_key=library_key, **c2l_kwargs,
            )
        else:
            proportions, cell_type_names, n_shared = _deconvolve_c2l(
                adata, ref_adata, batch_key_st=batch_key_st, **c2l_kwargs,
            )
        method_params = {
            "batch_key_ref": batch_key_ref,
            "batch_key_st": batch_key_st if not per_sample else None,
            "covariate_keys": covariate_keys,
            "N_cells_per_location": N_cells_per_location,
            "detection_alpha": detection_alpha,
            "max_epochs_ref": max_epochs_ref,
            "max_epochs_st": max_epochs_st,
            "batch_size_ref": batch_size_ref,
            "batch_size_st": batch_size_st,
            "num_samples_posterior": num_samples_posterior,
            "cell_count_cutoff": cell_count_cutoff,
            "cell_percentage_cutoff2": cell_percentage_cutoff2,
            "nonz_mean_cutoff": nonz_mean_cutoff,
        }
    else:
        raise ValueError(f"Unknown method={method!r}.")

    # ------------------------------------------------------------------ #
    # Write unified outputs to adata
    # ------------------------------------------------------------------ #
    adata.obsm["cell_type_abundances"] = proportions
    adata.obsm["q05_cell_abundance_w_sf"] = proportions   # back-compat alias

    for i, ct in enumerate(cell_type_names):
        adata.obs[ct] = proportions[:, i]

    dom_idx = np.argmax(proportions, axis=1)
    adata.obs["dominant_cell_type"] = pd.Categorical(
        [cell_type_names[i] for i in dom_idx],
        categories=cell_type_names,
    )

    # ------------------------------------------------------------------ #
    # Provenance
    # ------------------------------------------------------------------ #
    params = {
        "module": "spatial_deconvolve",
        "timestamp": datetime.now().isoformat(),
        "skipped": False,
        "params": {**common, **method_params},
        "outputs": {
            "method": method,
            "per_sample": per_sample,
            "library_key": library_key,
            "n_cell_types": len(cell_type_names),
            "cell_type_names": list(cell_type_names),
            "n_shared_genes": int(n_shared),
            "n_spots": int(adata.n_obs),
        },
    }
    adata.uns["omicsage_spatial_deconvolve"] = params
    return adata, params


# ---------------------------------------------------------------------------
# NNLS method
# ---------------------------------------------------------------------------


def _deconvolve_nnls(
    adata: ad.AnnData,
    ref_adata: ad.AnnData,
    cell_type_key: str,
    layer_ref: str,
    n_jobs: int,
    target_sum: float,
) -> tuple[np.ndarray, list, int]:
    """Non-negative least squares deconvolution.

    For each spot, solves ``argmin ||W p - x||^2  s.t. p >= 0`` where:
      - W is the reference signature matrix (genes × cell_types), built as
        the per-cell-type mean of library-size-normalised reference counts.
      - x is the spot's library-size-normalised expression vector.
      - p is the cell-type proportion vector (rescaled to sum to 1).

    Memory profile: O(n_genes × n_cell_types) for W plus O(n_genes) per spot.
    For Kuppe (~3,000 genes × ~14 cell types × ~11,700 spots) peak RAM is
    ~150 MB regardless of the spot count.
    """
    from scipy.optimize import nnls
    from joblib import Parallel, delayed

    # 1. Build per-cell-type reference signatures
    if layer_ref not in ref_adata.layers:
        raise ValueError(f"ref_adata.layers['{layer_ref}'] is missing.")

    cell_types = (
        ref_adata.obs[cell_type_key].astype("category").cat.categories.tolist()
    )
    ref_X = ref_adata.layers[layer_ref]
    ref_norm = _lib_normalize(ref_X, target_sum)

    sig = np.zeros((ref_adata.n_vars, len(cell_types)), dtype=np.float32)
    for i, ct in enumerate(cell_types):
        mask = (ref_adata.obs[cell_type_key].values == ct)
        if mask.sum() == 0:
            continue
        sig[:, i] = ref_norm[mask].mean(axis=0)

    ref_sig_df = pd.DataFrame(
        sig, index=ref_adata.var_names, columns=cell_types
    )

    # 2. Restrict to shared genes
    shared = [g for g in adata.var_names if g in ref_sig_df.index]
    if not shared:
        raise ValueError(
            "No shared genes between spatial and reference data. "
            "Check that both use the same gene ID format (e.g. ENSEMBL IDs)."
        )
    W = ref_sig_df.loc[shared].values.astype(np.float32)
    # Drop genes with zero signal in every cell type (no information)
    nonzero = W.sum(axis=1) > 0
    if not nonzero.all():
        W = W[nonzero]
        shared = [g for g, keep in zip(shared, nonzero) if keep]

    # 3. Library-size normalise spatial counts
    st_norm = _lib_normalize(
        adata[:, shared].layers["counts"], target_sum,
    ).astype(np.float32)

    # 4. NNLS per spot (parallel via threads — scipy.nnls releases the GIL)
    logger.info(
        "NNLS deconvolution: %d spots × %d cell types × %d shared genes",
        st_norm.shape[0], W.shape[1], W.shape[0],
    )

    def _solve(x):
        p, _ = nnls(W, x, maxiter=200)
        s = p.sum()
        return p / s if s > 0 else p

    proportions = np.asarray(
        Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_solve)(st_norm[i]) for i in range(st_norm.shape[0])
        ),
        dtype=np.float32,
    )
    return proportions, cell_types, len(shared)


def _lib_normalize(X, target_sum: float) -> np.ndarray:
    """Library-size normalise each row to ``target_sum``. Returns dense float32."""
    if sp.issparse(X):
        X = X.toarray().astype(np.float32)
    else:
        X = np.asarray(X, dtype=np.float32)
    row_sums = X.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    return (X / row_sums[:, None]) * target_sum


# ---------------------------------------------------------------------------
# cell2location method — split into reference-fit and spatial-fit steps
# ---------------------------------------------------------------------------


def _fit_c2l_reference(
    ref_adata: ad.AnnData,
    cell_type_key: str,
    layer_ref: str,
    batch_key_ref: Optional[str],
    covariate_keys: Optional[list],
    max_epochs_ref: int,
    batch_size_ref: int,
    cell_count_cutoff: int,
    cell_percentage_cutoff2: float,
    nonz_mean_cutoff: float,
    num_samples_posterior: int,
) -> tuple[pd.DataFrame, ad.AnnData]:
    """Fit cell2location RegressionModel — returns inferred signatures + filtered ref."""
    ref = ref_adata.copy()
    ref.X = ref.layers[layer_ref].copy()

    selected = c2l.utils.filtering.filter_genes(
        ref,
        cell_count_cutoff=cell_count_cutoff,
        cell_percentage_cutoff2=cell_percentage_cutoff2,
        nonz_mean_cutoff=nonz_mean_cutoff,
    )
    ref = ref[:, selected].copy()

    setup_kw: dict = dict(adata=ref, labels_key=cell_type_key, layer=layer_ref)
    if batch_key_ref and batch_key_ref in ref.obs.columns:
        setup_kw["batch_key"] = batch_key_ref
    if covariate_keys:
        valid = [k for k in covariate_keys if k in ref.obs.columns]
        if valid:
            setup_kw["categorical_covariate_keys"] = valid

    c2l.models.RegressionModel.setup_anndata(**setup_kw)
    reg_model = c2l.models.RegressionModel(ref)
    reg_model.train(
        max_epochs=max_epochs_ref,
        batch_size=batch_size_ref,
        train_size=1,
        lr=0.002,
    )
    reg_model.export_posterior(
        ref,
        sample_kwargs={
            "num_samples": num_samples_posterior,
            "batch_size": batch_size_ref,
        },
    )

    factor_names = ref.uns["mod"]["factor_names"]
    if "means_per_cluster_mu_fg" in ref.varm.keys():
        inf_aver = ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
        ].copy()
    else:
        inf_aver = ref.var[
            [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
        ].copy()
    inf_aver.columns = factor_names
    return inf_aver, ref


def _fit_c2l_spatial(
    adata_st: ad.AnnData,
    inf_aver: pd.DataFrame,
    batch_key_st: Optional[str],
    N_cells_per_location: int,
    detection_alpha: int,
    max_epochs_st: int,
    batch_size_st: Optional[int],
    num_samples_posterior: int,
) -> tuple[np.ndarray, list]:
    """Fit the cell2location spatial model and return (proportions, cell_type_names)."""
    setup_kw: dict = dict(adata=adata_st)
    if batch_key_st and batch_key_st in adata_st.obs.columns:
        setup_kw["batch_key"] = batch_key_st
    c2l.models.Cell2location.setup_anndata(**setup_kw)

    model = c2l.models.Cell2location(
        adata_st,
        cell_state_df=inf_aver,
        N_cells_per_location=N_cells_per_location,
        detection_alpha=detection_alpha,
    )
    model.train(
        max_epochs=max_epochs_st,
        batch_size=batch_size_st,
        train_size=1,
    )

    adata_st = model.export_posterior(
        adata_st,
        sample_kwargs={
            "num_samples": num_samples_posterior,
            "batch_size": batch_size_st or model.adata.n_obs,
        },
    )

    cell_type_names = list(adata_st.uns["mod"]["factor_names"])
    proportions = np.asarray(
        adata_st.obsm["q05_cell_abundance_w_sf"], dtype=np.float32
    )
    return proportions, cell_type_names


def _deconvolve_c2l(
    adata: ad.AnnData,
    ref_adata: ad.AnnData,
    cell_type_key: str,
    layer_ref: str,
    batch_key_ref: Optional[str],
    batch_key_st: Optional[str],
    covariate_keys: Optional[list],
    N_cells_per_location: int,
    detection_alpha: int,
    max_epochs_ref: int,
    max_epochs_st: int,
    batch_size_ref: int,
    batch_size_st: Optional[int],
    num_samples_posterior: int,
    cell_count_cutoff: int,
    cell_percentage_cutoff2: float,
    nonz_mean_cutoff: float,
) -> tuple[np.ndarray, list, int]:
    """Standard cell2location pipeline — all spots in one spatial fit."""
    inf_aver, ref = _fit_c2l_reference(
        ref_adata, cell_type_key, layer_ref,
        batch_key_ref, covariate_keys,
        max_epochs_ref, batch_size_ref,
        cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff,
        num_samples_posterior,
    )

    adata_st = adata.copy()
    adata_st.X = adata_st.layers["counts"].copy()
    shared = [g for g in adata_st.var_names if g in ref.var_names]
    if not shared:
        raise ValueError("No shared genes between spatial and reference data.")
    adata_st = adata_st[:, shared].copy()
    inf_aver_shared = inf_aver.loc[shared]

    proportions, cell_type_names = _fit_c2l_spatial(
        adata_st, inf_aver_shared,
        batch_key_st=batch_key_st,
        N_cells_per_location=N_cells_per_location,
        detection_alpha=detection_alpha,
        max_epochs_st=max_epochs_st,
        batch_size_st=batch_size_st,
        num_samples_posterior=num_samples_posterior,
    )
    return proportions, cell_type_names, len(shared)


def _deconvolve_c2l_per_sample(
    adata: ad.AnnData,
    ref_adata: ad.AnnData,
    library_key: str,
    cell_type_key: str,
    layer_ref: str,
    batch_key_ref: Optional[str],
    covariate_keys: Optional[list],
    N_cells_per_location: int,
    detection_alpha: int,
    max_epochs_ref: int,
    max_epochs_st: int,
    batch_size_ref: int,
    batch_size_st: Optional[int],
    num_samples_posterior: int,
    cell_count_cutoff: int,
    cell_percentage_cutoff2: float,
    nonz_mean_cutoff: float,
) -> tuple[np.ndarray, list, int]:
    """Per-sample cell2location — reference fit once, spatial fit per library.

    Each spatial fit holds only ~n_obs/N_libraries spots in memory, dropping
    peak RAM by roughly Nx compared to a single fit on all spots together.
    """
    inf_aver, ref = _fit_c2l_reference(
        ref_adata, cell_type_key, layer_ref,
        batch_key_ref, covariate_keys,
        max_epochs_ref, batch_size_ref,
        cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff,
        num_samples_posterior,
    )

    shared = [g for g in adata.var_names if g in ref.var_names]
    if not shared:
        raise ValueError("No shared genes between spatial and reference data.")
    inf_aver_shared = inf_aver.loc[shared]

    libraries = list(adata.obs[library_key].unique())
    logger.info(
        "Per-sample cell2location: %d libraries × spatial fits", len(libraries)
    )

    cell_type_names: Optional[list] = None
    proportions_all = np.zeros(
        (adata.n_obs, inf_aver_shared.shape[1]), dtype=np.float32
    )

    for i, lib_id in enumerate(libraries):
        mask = (adata.obs[library_key].values == lib_id)
        logger.info(
            "  [%d/%d] library=%s  spots=%d",
            i + 1, len(libraries), lib_id, int(mask.sum()),
        )
        adata_sub = adata[mask, shared].copy()
        adata_sub.X = adata_sub.layers["counts"].copy()
        # Inside a single library batch_key would have only one level — disable
        proportions, ct_names = _fit_c2l_spatial(
            adata_sub, inf_aver_shared,
            batch_key_st=None,
            N_cells_per_location=N_cells_per_location,
            detection_alpha=detection_alpha,
            max_epochs_st=max_epochs_st,
            batch_size_st=batch_size_st,
            num_samples_posterior=num_samples_posterior,
        )
        if cell_type_names is None:
            cell_type_names = ct_names
        proportions_all[mask] = proportions

    return proportions_all, cell_type_names, len(shared)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _skip_provenance(method: str, reason: str) -> dict:
    return {
        "module": "spatial_deconvolve",
        "timestamp": datetime.now().isoformat(),
        "skipped": True,
        "skip_reason": reason,
        "params": {"method": method},
        "outputs": {"method": method},
    }


def _check_c2l() -> None:
    if not _C2L_AVAILABLE:
        raise ImportError(
            "cell2location is required for method='cell2location'. "
            "Install with `pip install cell2location`, or use method='nnls' "
            "for a pure-Python alternative with no extra dependencies."
        )


def _validate_input(
    adata: ad.AnnData,
    ref_adata: Optional[ad.AnnData],
    cell_type_key: str,
    layer_ref: str,
    method: str,
) -> None:
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"Expected AnnData for adata, got {type(adata).__name__}")
    if adata.n_obs == 0:
        raise ValueError("adata has 0 observations.")
    if "counts" not in adata.layers:
        raise ValueError(
            "adata.layers['counts'] is missing. "
            "Run spatial_ingest() first."
        )
    if method not in ("nnls", "cell2location", "none"):
        raise ValueError(
            f"Unknown method={method!r}. "
            f"Expected 'nnls', 'cell2location', or 'none'."
        )
    if ref_adata is None or method == "none":
        return
    if not isinstance(ref_adata, ad.AnnData):
        raise TypeError(
            f"Expected AnnData for ref_adata, got {type(ref_adata).__name__}"
        )
    if ref_adata.n_obs == 0:
        raise ValueError("ref_adata has 0 observations.")
    if layer_ref not in ref_adata.layers:
        raise ValueError(
            f"ref_adata.layers['{layer_ref}'] is missing. "
            f"Available: {list(ref_adata.layers.keys())}"
        )
    if cell_type_key not in ref_adata.obs.columns:
        raise ValueError(
            f"cell_type_key='{cell_type_key}' not in ref_adata.obs. "
            f"Available: {list(ref_adata.obs.columns)}"
        )
