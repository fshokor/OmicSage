"""
cite_integration.py — Multi-modal CITE-seq integration for OmicSage.
 
Implements two fully-Python integration methods for paired RNA + ADT data:
 
  A. MOFA+  (Multi-Omics Factor Analysis)
     Linear factor model. Fast, interpretable, handles batch via groups_label.
     Reference: Argelaguet et al., Genome Biology 2020.
     Python: mu.tl.mofa(mdata, groups_label=batch_key)
 
  B. totalVI (Total Variational Inference)
     Deep generative model. Best scIB benchmark score on NeurIPS CITE-seq.
     Models RNA with NB distribution, ADT as NB foreground/background mixture.
     Requires raw counts for both RNA and ADT.
     Reference: Gayoso et al., Nature Methods 2021.
     Python: scvi.model.TOTALVI
 
  NOTE — WNN (Weighted Nearest Neighbor):
     WNN via muon.pp.neighbors is deferred to a future session.
     muon's internal Jaccard+Euclidean NNDescent computation hangs on small
     fixtures due to pynndescent behaviour.  The API contract is documented
     here for reference but the implementation is not included yet.
     WNN output keys (when implemented):
       mdata.obsm["X_umap_wnn"]          — WNN UMAP
       mdata.obsp["wnn_connectivities"]  — WNN graph connectivities
       mdata.obsp["wnn_distances"]       — WNN graph distances
 
Embedding keys written to mdata (OmicSage convention):
  mdata.obsm["X_mofa"]         — MOFA+ latent factors
  mdata.obsm["X_umap_mofa"]    — UMAP from MOFA+ embedding
  mdata.obsm["X_totalVI"]      — totalVI latent representation
  mdata.obsm["X_umap_totalVI"] — UMAP from totalVI embedding
 
Prerequisites per method
------------------------
MOFA+  : mdata["rna"].X  — log1p-normalized RNA
         mdata["adt"].X  — CLR-normalized ADT
         obs[batch_key] present on both modalities
totalVI: mdata["rna"].layers["counts"] — raw RNA integer counts
         mdata["adt"].layers["counts"] — raw ADT integer counts
         mdata["rna"].obs[batch_key]   — batch column
 
API
---
run_mofa(mdata, batch_key, n_factors=15, use_layer=None, random_state=0,
         inplace=False) -> (MuData, dict)
 
run_totalvi(mdata, batch_key, max_epochs=400, random_state=0,
            inplace=False) -> (MuData, dict)
 
References
----------
sc-best-practices ch.44:
  https://www.sc-best-practices.org/multimodal_integration/paired_integration.html
"""
 
from __future__ import annotations
 
import warnings
from datetime import datetime, timezone
from typing import Optional
 
import numpy as np
import scanpy as sc
from mudata import MuData
 
# Optional heavy dependencies — imported at module level so they are
# cached in sys.modules on first load and never re-imported mid-session.
# Each function still guards with a clear error if not installed.
try:
    import scvi as _scvi
    import scvi.model as _scvi_model
    import scvi.model._totalvi as _scvi_totalvi
    _TOTALVI = _scvi_totalvi.TOTALVI
    _SCVI_AVAILABLE = True
except (ModuleNotFoundError, ImportError, Exception):
    _SCVI_AVAILABLE = False
    _TOTALVI = None
 
 
# ---------------------------------------------------------------------------
# A. MOFA+
# ---------------------------------------------------------------------------
 
def run_mofa(
    mdata: MuData,
    batch_key: str,
    n_factors: int = 15,
    use_layer: Optional[str] = None,
    random_state: int = 0,
    inplace: bool = False,
) -> tuple[MuData, dict]:
    """
    Run MOFA+ integration on paired RNA + ADT data.
 
    Calls ``mu.tl.mofa`` with ``groups_label=batch_key`` for batch-aware
    factor decomposition.  The resulting latent factors are stored in
    ``mdata.obsm["X_mofa"]``.  A UMAP is computed from the MOFA+ embedding
    and stored in ``mdata.obsm["X_umap_mofa"]``.
 
    Parameters
    ----------
    mdata : MuData
        MuData with ``mdata["rna"]`` and ``mdata["adt"]`` modalities.
        - ``mdata["rna"].X``     — log1p-normalized RNA values
        - ``mdata["adt"].X``     — CLR-normalized ADT values
        - ``obs[batch_key]``     — batch column (present on both modalities)
    batch_key : str
        Column in ``.obs`` used as the group/batch variable for MOFA+.
        Typical value for GSE194122: ``"donor"``.
    n_factors : int
        Number of MOFA+ latent factors.  Default: 15.
    use_layer : str or None
        If set, uses this layer from each modality instead of ``.X``.
        Default: None (uses ``.X``).
    random_state : int
        Random seed passed to MOFA+.  Default: 0.
    inplace : bool
        If ``True``, modify ``mdata`` in place.
        If ``False`` (default), operate on a copy.
 
    Returns
    -------
    mdata_out : MuData
        MuData with:
          ``obsm["X_mofa"]``        — MOFA+ latent factors (n_cells × n_factors)
          ``obsm["X_umap_mofa"]``   — UMAP embedding (n_cells × 2)
          ``uns["omicsage_mofa"]``  — provenance dict
    metrics : dict
        Keys: n_cells, n_factors, batch_key, n_batches, batch_values,
              mofa_key, umap_key, umap_computed, random_state.
 
    Raises
    ------
    ImportError
        If ``muon`` or ``mofapy2`` are not installed.
    KeyError
        If required modalities or batch_key are missing.
    """
    try:
        import muon
    except ImportError as exc:
        raise ImportError(
            "muon is required for MOFA+ integration. "
            "Install with: pip install muon"
        ) from exc
 
    try:
        import mofapy2  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "mofapy2 is required for MOFA+ integration. "
            "Install with: pip install mofapy2"
        ) from exc
 
    _validate_mofa_inputs(mdata, batch_key)
 
    mdata_out = mdata if inplace else mdata.copy()
 
    # muon namespaces conflicting obs columns as "rna:donor" / "adt:donor"
    # in mdata.obs.  MOFA+ reads groups_label from mdata.obs directly, so
    # we push the batch column to the top-level mdata.obs before calling
    # mu.tl.mofa.  Leaving it there is intentional — useful for downstream
    # visualisations.
    mdata_out.obs[batch_key] = mdata_out["rna"].obs[batch_key].values
 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        muon.tl.mofa(
            mdata_out,
            groups_label=batch_key,
            use_layer=use_layer,
            n_factors=n_factors,
            seed=random_state,
            quiet=True,
        )
 
    # Compute UMAP from MOFA+ embedding and rename to avoid collisions
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(mdata_out, use_rep="X_mofa", random_state=random_state)
        sc.tl.umap(mdata_out, random_state=random_state)
 
    mdata_out.obsm["X_umap_mofa"] = mdata_out.obsm.pop("X_umap")
 
    batch_values = sorted(
        mdata_out["rna"].obs[batch_key].astype(str).unique().tolist()
    )
    n_factors_actual = mdata_out.obsm["X_mofa"].shape[1]
 
    metrics: dict = {
        "n_cells": int(mdata_out.n_obs),
        "n_factors": n_factors_actual,
        "batch_key": batch_key,
        "n_batches": len(batch_values),
        "batch_values": batch_values,
        "mofa_key": "X_mofa",
        "umap_key": "X_umap_mofa",
        "umap_computed": True,
        "random_state": random_state,
    }
 
    mdata_out.uns["omicsage_mofa"] = {
        "module": "cite_integration.run_mofa",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "batch_key": batch_key,
            "n_factors": n_factors,
            "use_layer": use_layer,
            "random_state": random_state,
        },
        "outputs": {
            "mofa_key": "X_mofa",
            "umap_key": "X_umap_mofa",
            "n_factors_actual": n_factors_actual,
        },
    }
 
    return mdata_out, metrics
 
 
# ---------------------------------------------------------------------------
# B. totalVI
# ---------------------------------------------------------------------------
 
def run_totalvi(
    mdata: MuData,
    batch_key: str,
    max_epochs: int = 400,
    random_state: int = 0,
    inplace: bool = False,
) -> tuple[MuData, dict]:
    """
    Run totalVI integration on paired RNA + ADT data.
 
    Constructs a single AnnData with RNA counts in ``.X`` and ADT counts in
    ``obsm["protein_expression"]``, as required by ``scvi.model.TOTALVI``.
    The latent representation is stored in ``mdata.obsm["X_totalVI"]``.
    A UMAP is computed from the totalVI embedding and stored in
    ``mdata.obsm["X_umap_totalVI"]``.
 
    totalVI models RNA counts with a negative-binomial distribution and ADT
    counts as a foreground/background NB mixture, making it the only method
    here that explicitly accounts for ADT background noise.
 
    Parameters
    ----------
    mdata : MuData
        MuData with ``mdata["rna"]`` and ``mdata["adt"]`` modalities.
        - ``mdata["rna"].layers["counts"]`` — raw RNA integer counts (required)
        - ``mdata["adt"].layers["counts"]`` — raw ADT integer counts (required)
        - ``mdata["rna"].obs[batch_key]``   — batch column on RNA modality
    batch_key : str
        Column in ``mdata["rna"].obs`` used as batch variable.
        Typical value for GSE194122: ``"donor"``.
    max_epochs : int
        Maximum training epochs.  Default: 400.
        Use a smaller value (e.g. 2) for fast tests.
    random_state : int
        Random seed for scvi-tools.  Default: 0.
    inplace : bool
        If ``True``, modify ``mdata`` in place.
        If ``False`` (default), operate on a copy.
 
    Returns
    -------
    mdata_out : MuData
        MuData with:
          ``obsm["X_totalVI"]``       — totalVI latent representation
                                        (n_cells × n_latent)
          ``obsm["X_umap_totalVI"]``  — UMAP embedding (n_cells × 2)
          ``uns["omicsage_totalVI"]`` — provenance dict
    metrics : dict
        Keys: n_cells, n_latent, batch_key, n_batches, batch_values,
              max_epochs, totalvi_key, umap_key, umap_computed, random_state.
 
    Raises
    ------
    ImportError
        If ``scvi-tools`` is not installed.
    KeyError
        If required modalities, layers, or batch_key are missing.
    """
    _validate_totalvi_inputs(mdata, batch_key)

    try:
        import scvi
        from scvi.model._totalvi import TOTALVI
    except (ModuleNotFoundError, ImportError) as exc:
        raise ImportError(
            "scvi-tools is required for totalVI integration. "
            "Install with: pip install scvi-tools"
        ) from exc
 
    mdata_out = mdata if inplace else mdata.copy()
 
    rna = mdata_out["rna"]
    adt = mdata_out["adt"]
 
    import anndata as ad
    import scipy.sparse as sp
 
    rna_counts = rna.layers["counts"]
    adt_counts = adt.layers["counts"]
 
    # totalVI expects protein_expression as a dense array or DataFrame
    if sp.issparse(adt_counts):
        adt_dense = adt_counts.toarray().astype(np.float32)
    else:
        adt_dense = np.asarray(adt_counts, dtype=np.float32)
 
    adata_vi = ad.AnnData(X=rna_counts.copy())
    adata_vi.obs_names = rna.obs_names.copy()
    adata_vi.var_names = rna.var_names.copy()
    adata_vi.obs[batch_key] = rna.obs[batch_key].copy()
    adata_vi.obsm["protein_expression"] = adt_dense
    adata_vi.layers["counts"] = rna_counts.copy()
 
    # scvi.settings available from scvi-tools >= 0.16; guard for older envs
    if hasattr(scvi, "settings"):
        scvi.settings.seed = random_state
 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        TOTALVI.setup_anndata(
            adata_vi,
            protein_expression_obsm_key="protein_expression",
            layer="counts",
            batch_key=batch_key,
        )
        vae = TOTALVI(adata_vi)
        vae.train(max_epochs=max_epochs)
 
    latent = vae.get_latent_representation()
    mdata_out.obsm["X_totalVI"] = latent
 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(mdata_out, use_rep="X_totalVI", random_state=random_state)
        sc.tl.umap(mdata_out, random_state=random_state)
 
    mdata_out.obsm["X_umap_totalVI"] = mdata_out.obsm.pop("X_umap")
 
    batch_values = sorted(rna.obs[batch_key].astype(str).unique().tolist())
    n_latent = latent.shape[1]
 
    metrics: dict = {
        "n_cells": int(mdata_out.n_obs),
        "n_latent": n_latent,
        "batch_key": batch_key,
        "n_batches": len(batch_values),
        "batch_values": batch_values,
        "max_epochs": max_epochs,
        "totalvi_key": "X_totalVI",
        "umap_key": "X_umap_totalVI",
        "umap_computed": True,
        "random_state": random_state,
    }
 
    mdata_out.uns["omicsage_totalVI"] = {
        "module": "cite_integration.run_totalvi",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "batch_key": batch_key,
            "max_epochs": max_epochs,
            "random_state": random_state,
        },
        "outputs": {
            "totalvi_key": "X_totalVI",
            "umap_key": "X_umap_totalVI",
            "n_latent": n_latent,
        },
    }
 
    return mdata_out, metrics
 
 
# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
 
def _validate_mofa_inputs(mdata: MuData, batch_key: str) -> None:
    for mod in ("rna", "adt"):
        if mod not in mdata.mod:
            raise KeyError(
                f"mdata must contain a '{mod}' modality. "
                f"Found: {list(mdata.mod.keys())}"
            )
    if batch_key not in mdata["rna"].obs.columns:
        raise KeyError(
            f"batch_key='{batch_key}' not found in mdata['rna'].obs. "
            f"Available: {list(mdata['rna'].obs.columns)}"
        )
 
 
def _validate_totalvi_inputs(mdata: MuData, batch_key: str) -> None:
    for mod in ("rna", "adt"):
        if mod not in mdata.mod:
            raise KeyError(
                f"mdata must contain a '{mod}' modality. "
                f"Found: {list(mdata.mod.keys())}"
            )
    if "counts" not in mdata["rna"].layers:
        raise KeyError(
            "mdata['rna'].layers['counts'] not found. "
            "totalVI requires raw RNA counts. "
            "Ensure ingest.py preserved raw counts in layers['counts']."
        )
    if "counts" not in mdata["adt"].layers:
        raise KeyError(
            "mdata['adt'].layers['counts'] not found. "
            "totalVI requires raw ADT counts. "
            "Ensure adt_normalize.py preserved raw counts in layers['counts']."
        )
    if batch_key not in mdata["rna"].obs.columns:
        raise KeyError(
            f"batch_key='{batch_key}' not found in mdata['rna'].obs. "
            f"Available: {list(mdata['rna'].obs.columns)}"
        )