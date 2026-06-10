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
    log_every_n_epochs: int = 10,
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
            log_every_n_epochs=log_every_n_epochs,
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

    For each spot solves ``argmin ||W p - x||^2  s.t. p >= 0`` where:
      - W  is the (n_shared_genes × n_cell_types) signature matrix
      - x  is the spot's library-size-normalised expression vector
      - p  is the cell-type proportion vector (rescaled to sum to 1)

    Memory design
    -------------
    The reference can be large (e.g. 50k cells × 36k genes = 7 GB dense).
    We avoid ever densifying it by:
      1. Finding shared genes first (index lookup, no data loaded)
      2. Subsetting the sparse reference to shared genes only (~3k)
      3. Computing per-cell-type means one cell type at a time from the
         sparse subset — peak allocation is
         max(n_cells_per_celltype) × n_shared × 4 B ≈ 40–80 MB
      4. The spatial matrix is small (11k × 3k × 4 B ≈ 130 MB) and is
         the only dense allocation

    Combined peak RAM for Kuppe: ~200–300 MB.
    """
    from scipy.optimize import nnls
    from joblib import Parallel, delayed

    if layer_ref not in ref_adata.layers:
        raise ValueError(f"ref_adata.layers['{layer_ref}'] is missing.")

    cell_types = (
        ref_adata.obs[cell_type_key].astype("category").cat.categories.tolist()
    )

    # ── Step 1: find shared genes while everything is still sparse ──────────
    spatial_genes = set(adata.var_names)
    shared = [g for g in ref_adata.var_names if g in spatial_genes]
    if not shared:
        raise ValueError(
            "No shared genes between spatial and reference data. "
            "Check that both use the same gene ID format (e.g. ENSEMBL IDs)."
        )
    logger.info(
        "NNLS: %d shared genes, %d cell types, %d ref cells, %d spots",
        len(shared), len(cell_types), ref_adata.n_obs, adata.n_obs,
    )

    # ── Step 2: build signature matrix — sparse, one cell type at a time ────
    # ref_sub is still SPARSE: n_ref_cells × n_shared — never fully densified
    shared_idx = {g: i for i, g in enumerate(ref_adata.var_names)}
    shared_pos = [shared_idx[g] for g in shared]

    sig = np.zeros((len(shared), len(cell_types)), dtype=np.float32)
    ct_labels = np.asarray(ref_adata.obs[cell_type_key].values)

    for i, ct in enumerate(cell_types):
        mask = (ct_labels == ct)
        n_ct = mask.sum()
        if n_ct == 0:
            continue
        # Slice sparse: n_ct_cells × n_shared — max ~3.5k × 3k × 4B ≈ 42 MB
        ct_X = ref_adata.layers[layer_ref][mask]
        if sp.issparse(ct_X):
            ct_X = ct_X[:, shared_pos]
        else:
            ct_X = np.asarray(ct_X, dtype=np.float32)[:, shared_pos]
        # Per-row library-size normalisation (sparse-aware)
        row_sums = np.array(
            ct_X.sum(axis=1) if sp.issparse(ct_X)
            else ct_X.sum(axis=1)
        ).ravel().astype(np.float32)
        row_sums[row_sums == 0] = 1.0
        if sp.issparse(ct_X):
            ct_norm = ct_X.multiply(
                (target_sum / row_sums)[:, None]
            ).toarray().astype(np.float32)
        else:
            ct_norm = (ct_X / row_sums[:, None]) * target_sum
        sig[:, i] = ct_norm.mean(axis=0)
        del ct_X, ct_norm  # free immediately

    # Drop genes with zero signal across all cell types
    nonzero = sig.sum(axis=1) > 0
    if not nonzero.all():
        sig    = sig[nonzero]
        shared = [g for g, keep in zip(shared, nonzero) if keep]
    W = sig  # (n_shared, n_cell_types) — tiny

    # ── Step 3: normalise spatial counts (small matrix) ─────────────────────
    # adata[:, shared] subset: 11k × 3k, dense float32 ≈ 130 MB
    spatial_pos = [list(adata.var_names).index(g) for g in shared]
    st_counts = adata.layers["counts"]
    if sp.issparse(st_counts):
        st_sub = st_counts[:, spatial_pos].toarray().astype(np.float32)
    else:
        st_sub = np.asarray(st_counts, dtype=np.float32)[:, spatial_pos]
    row_sums = st_sub.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    st_norm = (st_sub / row_sums[:, None]) * target_sum
    del st_sub  # free raw counts copy

    # ── Step 4: NNLS per spot (threads share W and st_norm, no copies) ──────
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


def _make_epoch_callback(label: str, max_epochs: int, log_every_n: int):
    """Lightning callback: log epoch/loss/ETA every log_every_n epochs.

    Tries lightning.pytorch, falls back to pytorch_lightning.
    Returns None if neither is available.
    """
    import time as _time
    try:
        try:
            from lightning.pytorch.callbacks import Callback
        except ImportError:
            from pytorch_lightning.callbacks import Callback
    except ImportError:
        return None

    _nd = len(str(max_epochs))

    class _EpochLog(Callback):
        def on_train_start(self, trainer, pl_module):
            self._t0 = _time.time()
            print(f"    [{label}] starting -- {max_epochs} epochs", flush=True)

        def on_train_epoch_end(self, trainer, pl_module):
            epoch = trainer.current_epoch + 1
            if epoch % log_every_n != 0 and epoch != max_epochs:
                return
            elapsed = _time.time() - self._t0
            loss_val = trainer.callback_metrics.get("train_loss_epoch")
            loss_str = (
                f"  loss={float(loss_val):.3f}" if loss_val is not None else ""
            )
            eta_str = ""
            if epoch > 0:
                eta_s = ((max_epochs - epoch) * elapsed) / epoch
                eta_str = (
                    f"  ETA {eta_s / 60:.1f}min" if eta_s > 5 else "  ETA <1min"
                )
            print(
                f"    [{label}] epoch {epoch:{_nd}d}/{max_epochs}{loss_str}{eta_str}",
                flush=True,
            )

        def on_train_end(self, trainer, pl_module):
            elapsed = _time.time() - self._t0
            print(f"    [{label}] done -- {elapsed / 60:.1f} min", flush=True)

    return _EpochLog()


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
    log_every_n_epochs: int = 10,
) -> tuple[pd.DataFrame, ad.AnnData]:
    """Fit cell2location RegressionModel — returns inferred signatures + filtered ref."""
    ref = ref_adata.copy()
    ref.X = ref.layers[layer_ref].copy()

    # filter_genes internally calls log10(nonz_mean) which raises divide-by-zero
    # for genes with zero mean — those genes are discarded by the filter anyway,
    # so the warning is expected and harmless. Suppress at both numpy and Python
    # warning levels because pandas re-raises numpy RuntimeWarnings via its own
    # arraylike dispatch path.
    import warnings as _warnings
    with _warnings.catch_warnings(), np.errstate(divide="ignore", invalid="ignore"):
        _warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")
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
    _cb = _make_epoch_callback("ref model", max_epochs_ref, log_every_n_epochs)
    _kw: dict = dict(
        max_epochs=max_epochs_ref, batch_size=batch_size_ref,
        train_size=1, lr=0.002, enable_progress_bar=False,
    )
    if _cb is not None:
        _kw["callbacks"] = [_cb]
    reg_model.train(**_kw)
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
    label: str = "spatial",
    log_every_n_epochs: int = 10,
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
    _cb = _make_epoch_callback(label, max_epochs_st, log_every_n_epochs)
    _kw: dict = dict(
        max_epochs=max_epochs_st, batch_size=batch_size_st,
        train_size=1, enable_progress_bar=False,
    )
    if _cb is not None:
        _kw["callbacks"] = [_cb]
    model.train(**_kw)

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
    import time as _time
    t0 = _time.time()
    print(
        f"  [deconvolve] step 1/2 -- reference model ({max_epochs_ref} epochs) ...",
        flush=True,
    )
    inf_aver, ref = _fit_c2l_reference(
        ref_adata, cell_type_key, layer_ref,
        batch_key_ref, covariate_keys,
        max_epochs_ref, batch_size_ref,
        cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff,
        num_samples_posterior,
        log_every_n_epochs=log_every_n_epochs,
    )
    print(
        f"  [deconvolve] reference done in {(_time.time() - t0) / 60:.1f} min",
        flush=True,
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
    log_every_n_epochs: int = 10,
) -> tuple[np.ndarray, list, int]:
    """Per-sample cell2location — reference fit once, spatial fit per library.

    Each spatial fit holds only ~n_obs/N_libraries spots in memory, dropping
    peak RAM by roughly Nx compared to a single fit on all spots together.
    """
    import time as _time
    t0 = _time.time()
    print(
        f"  [deconvolve] step 1/2 -- reference model ({max_epochs_ref} epochs) ...",
        flush=True,
    )
    inf_aver, ref = _fit_c2l_reference(
        ref_adata, cell_type_key, layer_ref,
        batch_key_ref, covariate_keys,
        max_epochs_ref, batch_size_ref,
        cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff,
        num_samples_posterior,
        log_every_n_epochs=log_every_n_epochs,
    )
    print(
        f"  [deconvolve] reference done in {(_time.time() - t0) / 60:.1f} min",
        flush=True,
    )

    shared = [g for g in adata.var_names if g in ref.var_names]
    if not shared:
        raise ValueError("No shared genes between spatial and reference data.")
    inf_aver_shared = inf_aver.loc[shared]

    libraries = list(adata.obs[library_key].unique())
    print(
        f"  [deconvolve] step 2/2 -- spatial fits: "
        f"{len(libraries)} libraries x {max_epochs_st} epochs each",
        flush=True,
    )

    cell_type_names: Optional[list] = None
    proportions_all = np.zeros(
        (adata.n_obs, inf_aver_shared.shape[1]), dtype=np.float32
    )

    for i, lib_id in enumerate(libraries):
        mask = (adata.obs[library_key].values == lib_id)
        n_spots = int(mask.sum())
        t_lib = _time.time()
        print(
            f"  [deconvolve]   [{i + 1}/{len(libraries)}] {lib_id} ({n_spots:,} spots) ...",
            flush=True,
        )
        adata_sub = adata[mask, shared].copy()
        adata_sub.X = adata_sub.layers["counts"].copy()
        proportions, ct_names = _fit_c2l_spatial(
            adata_sub, inf_aver_shared,
            batch_key_st=None,
            N_cells_per_location=N_cells_per_location,
            detection_alpha=detection_alpha,
            max_epochs_st=max_epochs_st,
            batch_size_st=batch_size_st,
            num_samples_posterior=num_samples_posterior,
            label=f"{lib_id} [{i + 1}/{len(libraries)}]",
            log_every_n_epochs=log_every_n_epochs,
        )
        print(
            f"  [deconvolve]   [{i + 1}/{len(libraries)}] {lib_id} done "
            f"in {(_time.time() - t_lib) / 60:.1f} min",
            flush=True,
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
