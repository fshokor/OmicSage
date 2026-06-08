"""
spatial_deconvolve.py — OmicSage Phase 7, Session 4
Cell type deconvolution of Visium spots using cell2location.

Input:  AnnData from spatial_cluster (has layers["counts"], obsm["spatial"])
        Optional ref_adata: paired scRNA-seq reference AnnData
Output: AnnData with obsm["q05_cell_abundance_w_sf"] (cell type abundances),
        obs columns per cell type (5% quantile abundances),
        provenance in uns["omicsage_spatial_deconvolve"]

Two modes:
  ref_adata provided → run full cell2location pipeline
  ref_adata is None  → skip, write skipped=True in provenance

Pipeline (when ref_adata provided):
  1. Prepare reference: reset X to raw counts, gene selection
  2. Fit RegressionModel on reference → cell type signatures (inf_aver)
  3. Subset spatial + reference to shared genes
  4. Fit Cell2location on spatial data → cell type abundances per spot
  5. Export posterior → obsm["q05_cell_abundance_w_sf"]
  6. Write per-cell-type columns to obs (5% quantile)

Notes
-----
- All cell2location hyperparameters are arguments — different tissues
  require different N_cells_per_location, batch keys, etc.
- cell_type_key, batch_key_ref, batch_key_st, covariate_keys are all
  configurable to support any paired dataset.
- Device (GPU/CPU) is detected automatically by scvi-tools/PyTorch.
  No explicit use_gpu flag needed — falls back to CPU if no GPU present.
- Gene ID format: spatial var_names and reference var_names must match.
  For Kuppe data this means both use ENSEMBL IDs after h5ad ingest swap.
"""

from __future__ import annotations

from datetime import datetime
from typing import Optional

import anndata as ad
import numpy as np
import scanpy as sc

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
    cell_type_key: str = "cell_type_original",
    batch_key_ref: Optional[str] = "donor_id",
    batch_key_st: Optional[str] = "patient",
    covariate_keys: Optional[list] = None,
    layer_ref: str = "counts",
    N_cells_per_location: int = 8,
    detection_alpha: int = 20,
    max_epochs_ref: int = 250,
    max_epochs_st: int = 30000,
    batch_size_ref: int = 2500,
    cell_count_cutoff: int = 5,
    cell_percentage_cutoff2: float = 0.03,
    nonz_mean_cutoff: float = 1.12,
    inplace: bool = False,
) -> tuple[ad.AnnData, dict]:
    """Deconvolve Visium spots into cell type abundances using cell2location.

    Parameters
    ----------
    adata
        AnnData produced by :func:`spatial_cluster`.
        Must contain ``layers["counts"]`` with raw integer counts.
    ref_adata
        Paired scRNA-seq reference AnnData.
        When ``None``, deconvolution is skipped and provenance records
        ``skipped=True`` with a reason message.
    cell_type_key
        Column in ``ref_adata.obs`` containing cell type labels.
        Default: ``"cell_type_original"`` (Kuppe dataset).
    batch_key_ref
        Batch key in ``ref_adata.obs`` for the reference model.
        ``None`` disables batch correction in the reference model.
        Default: ``"donor_id"`` (Kuppe dataset).
    batch_key_st
        Batch key in ``adata.obs`` for the spatial model.
        ``None`` disables batch correction in the spatial model.
        Default: ``"patient"`` (Kuppe dataset).
    covariate_keys
        List of categorical covariate keys in ``ref_adata.obs``
        to include in the reference model.
        Default: ``None`` (no covariates).
    layer_ref
        Layer in ``ref_adata`` containing raw counts.
        Default: ``"counts"``.
    N_cells_per_location
        Expected number of cells per Visium spot.
        Obtain from paired histology images.
        Default: ``8`` (appropriate for heart tissue).
    detection_alpha
        Technical variability in RNA detection sensitivity.
        Higher = less regularisation. Default: ``20``.
    max_epochs_ref
        Training epochs for the reference RegressionModel.
        Default: ``250``.
    max_epochs_st
        Training epochs for the Cell2location spatial model.
        Default: ``30000``.
    batch_size_ref
        Batch size for reference model training.
        Default: ``2500``.
    cell_count_cutoff
        Minimum number of cells expressing a gene (gene selection).
        Default: ``5``.
    cell_percentage_cutoff2
        Minimum percentage of cells expressing a gene (gene selection).
        Default: ``0.03``.
    nonz_mean_cutoff
        Minimum mean expression in expressing cells (gene selection).
        Default: ``1.12``.
    inplace
        If ``False`` (default), operate on a copy.
        If ``True``, modify *adata* in place.

    Returns
    -------
    adata
        AnnData with (when ref_adata provided):

        - ``obsm["q05_cell_abundance_w_sf"]`` — cell type abundances (n_obs × n_cell_types)
        - ``obs[cell_type]`` columns — per-cell-type 5% quantile abundance
        - ``uns["omicsage_spatial_deconvolve"]`` — provenance

    params
        Provenance dictionary.

    Raises
    ------
    ImportError
        If cell2location is not installed and ref_adata is provided.
    TypeError
        If adata is not an AnnData object.
    ValueError
        If required keys are missing from adata or ref_adata.
    """
    _validate_input(adata, ref_adata, cell_type_key, layer_ref)

    adata = adata if inplace else adata.copy()

    # ------------------------------------------------------------------ #
    # Skip gracefully when no reference provided
    # ------------------------------------------------------------------ #
    if ref_adata is None:
        params = {
            "module": "spatial_deconvolve",
            "timestamp": datetime.now().isoformat(),
            "skipped": True,
            "skip_reason": (
                "ref_adata=None — deconvolution requires a paired scRNA-seq "
                "reference. Provide ref_adata to run cell2location."
            ),
            "params": {},
            "outputs": {},
        }
        adata.uns["omicsage_spatial_deconvolve"] = params
        return adata, params

    _check_c2l()

    # ------------------------------------------------------------------ #
    # 1. Prepare reference: set X to raw counts
    # ------------------------------------------------------------------ #
    ref = ref_adata.copy()
    ref.X = ref.layers[layer_ref].copy()

    # ------------------------------------------------------------------ #
    # 2. Gene selection on reference (permissive thresholds)
    # ------------------------------------------------------------------ #
    selected = c2l.utils.filtering.filter_genes(
        ref,
        cell_count_cutoff=cell_count_cutoff,
        cell_percentage_cutoff2=cell_percentage_cutoff2,
        nonz_mean_cutoff=nonz_mean_cutoff,
    )
    ref = ref[:, selected].copy()

    # ------------------------------------------------------------------ #
    # 3. Subset spatial data to shared genes
    # ------------------------------------------------------------------ #
    # Spatial data uses raw counts from layers["counts"]
    adata_st = adata.copy()
    adata_st.X = adata_st.layers["counts"].copy()

    shared_genes = [g for g in adata_st.var_names if g in ref.var_names]
    if len(shared_genes) == 0:
        raise ValueError(
            "No shared genes between spatial and reference data. "
            "Check that both use the same gene ID format (e.g. ENSEMBL IDs)."
        )
    ref = ref[:, shared_genes].copy()
    adata_st = adata_st[:, shared_genes].copy()

    # ------------------------------------------------------------------ #
    # 4. Fit RegressionModel on reference → cell type signatures
    # ------------------------------------------------------------------ #
    setup_ref_kwargs: dict = dict(
        adata=ref,
        labels_key=cell_type_key,
        layer=layer_ref,
    )
    if batch_key_ref is not None and batch_key_ref in ref.obs.columns:
        setup_ref_kwargs["batch_key"] = batch_key_ref
    if covariate_keys:
        valid_covs = [k for k in covariate_keys if k in ref.obs.columns]
        if valid_covs:
            setup_ref_kwargs["categorical_covariate_keys"] = valid_covs

    c2l.models.RegressionModel.setup_anndata(**setup_ref_kwargs)

    reg_model = c2l.models.RegressionModel(ref)
    reg_model.train(
        max_epochs=max_epochs_ref,
        batch_size=batch_size_ref,
        train_size=1,
        lr=0.002,
    )

    reg_model.export_posterior(
        ref,
        sample_kwargs={"num_samples": 1000, "batch_size": batch_size_ref},
    )

    # Extract inferred cell type signatures
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

    # ------------------------------------------------------------------ #
    # 5. Fit Cell2location on spatial data
    # ------------------------------------------------------------------ #
    setup_st_kwargs: dict = dict(adata=adata_st)
    if batch_key_st is not None and batch_key_st in adata_st.obs.columns:
        setup_st_kwargs["batch_key"] = batch_key_st

    c2l.models.Cell2location.setup_anndata(**setup_st_kwargs)

    c2l_model = c2l.models.Cell2location(
        adata_st,
        cell_state_df=inf_aver,
        N_cells_per_location=N_cells_per_location,
        detection_alpha=detection_alpha,
    )
    c2l_model.train(
        max_epochs=max_epochs_st,
        batch_size=None,
        train_size=1,
    )

    # ------------------------------------------------------------------ #
    # 6. Export posterior → cell type abundances
    # ------------------------------------------------------------------ #
    adata_st = c2l_model.export_posterior(
        adata_st,
        sample_kwargs={
            "num_samples": 1000,
            "batch_size": c2l_model.adata.n_obs,
        },
    )

    # Write 5% quantile abundances to obs columns (one per cell type)
    cell_type_names = adata_st.uns["mod"]["factor_names"]
    adata_st.obs[cell_type_names] = adata_st.obsm["q05_cell_abundance_w_sf"]

    # Copy deconvolution results back to the main adata (shared obs index)
    adata.obsm["q05_cell_abundance_w_sf"] = adata_st.obsm[
        "q05_cell_abundance_w_sf"
    ]
    for ct in cell_type_names:
        adata.obs[ct] = adata_st.obs[ct]
    if "mod" in adata_st.uns:
        adata.uns["cell2location_mod"] = adata_st.uns["mod"]

    # ------------------------------------------------------------------ #
    # 7. Provenance
    # ------------------------------------------------------------------ #
    n_shared = len(shared_genes)
    params = {
        "module": "spatial_deconvolve",
        "timestamp": datetime.now().isoformat(),
        "skipped": False,
        "params": {
            "cell_type_key": cell_type_key,
            "batch_key_ref": batch_key_ref,
            "batch_key_st": batch_key_st,
            "covariate_keys": covariate_keys,
            "layer_ref": layer_ref,
            "N_cells_per_location": N_cells_per_location,
            "detection_alpha": detection_alpha,
            "max_epochs_ref": max_epochs_ref,
            "max_epochs_st": max_epochs_st,
            "batch_size_ref": batch_size_ref,
            "cell_count_cutoff": cell_count_cutoff,
            "cell_percentage_cutoff2": cell_percentage_cutoff2,
            "nonz_mean_cutoff": nonz_mean_cutoff,
        },
        "outputs": {
            "n_cell_types": len(cell_type_names),
            "cell_type_names": list(cell_type_names),
            "n_shared_genes": n_shared,
            "n_genes_after_selection": int(selected.sum()),
            "n_spots": int(adata.n_obs),
        },
    }
    adata.uns["omicsage_spatial_deconvolve"] = params
    return adata, params


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _check_c2l() -> None:
    if not _C2L_AVAILABLE:
        raise ImportError(
            "cell2location is required for spatial deconvolution. "
            "Install with: pip install cell2location"
        )


def _validate_input(
    adata: ad.AnnData,
    ref_adata: Optional[ad.AnnData],
    cell_type_key: str,
    layer_ref: str,
) -> None:
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"Expected AnnData for adata, got {type(adata).__name__}")
    if adata.n_obs == 0:
        raise ValueError("adata has 0 observations. Cannot run deconvolution.")
    if "counts" not in adata.layers:
        raise ValueError(
            "adata.layers['counts'] is missing. "
            "Run spatial_ingest() and spatial_reduce() first."
        )
    if ref_adata is None:
        return  # graceful skip — no further validation needed
    if not isinstance(ref_adata, ad.AnnData):
        raise TypeError(
            f"Expected AnnData for ref_adata, got {type(ref_adata).__name__}"
        )
    if ref_adata.n_obs == 0:
        raise ValueError("ref_adata has 0 observations.")
    if layer_ref not in ref_adata.layers:
        raise ValueError(
            f"ref_adata.layers['{layer_ref}'] is missing. "
            f"Available layers: {list(ref_adata.layers.keys())}"
        )
    if cell_type_key not in ref_adata.obs.columns:
        raise ValueError(
            f"cell_type_key='{cell_type_key}' not found in ref_adata.obs. "
            f"Available columns: {list(ref_adata.obs.columns)}"
        )
