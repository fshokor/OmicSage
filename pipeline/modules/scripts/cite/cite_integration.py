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

  C. Both (run_both)
     Runs MOFA+ then totalVI sequentially on the same MuData.
     Both embeddings and UMAPs are written; scib metrics compare them.

scib benchmark metrics
----------------------
When ``compute_scib=True`` is passed to any run_* function, scib metrics are
computed for the resulting embedding and stored in ``metrics["scib"]``.

Metrics computed (matching sc-best-practices paired integration tutorial):
  Batch correction:
    iLISI      — higher is better (graph-based batch mixing)
    graph_conn — higher is better (global connectivity)
    kBET       — higher is better (local batch correction; skipped if slow)
  Bio conservation:
    cLISI      — higher is better (cell-type mixing)
    ASW_label  — higher is better (silhouette for cell types)

scib package: pip install scib

NOTE — WNN (Weighted Nearest Neighbor):
   WNN via muon.pp.neighbors is deferred to a future session.

Embedding keys written to mdata (OmicSage convention):
  mdata.obsm["X_mofa"]         — MOFA+ latent factors
  mdata.obsm["X_umap_mofa"]    — UMAP from MOFA+ embedding
  mdata.obsm["X_totalVI"]      — totalVI latent representation
  mdata.obsm["X_umap_totalVI"] — UMAP from totalVI embedding

API
---
run_mofa(mdata, batch_key, n_factors=15, use_layer=None, random_state=0,
         inplace=False, compute_scib=False, cell_type_key=None) -> (MuData, dict)

run_totalvi(mdata, batch_key, max_epochs=400, random_state=0,
            inplace=False, compute_scib=False, cell_type_key=None) -> (MuData, dict)

run_both(mdata, batch_key, n_factors=15, max_epochs=400, random_state=0,
         inplace=False, compute_scib=False, cell_type_key=None) -> (MuData, dict)

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
# scib metrics helper
# ---------------------------------------------------------------------------

def _compute_scib_metrics(
    mdata: MuData,
    embed_key: str,
    batch_key: str,
    cell_type_key: Optional[str],
    label: str = "",
) -> dict:
    """
    Compute scib integration benchmark metrics for a given embedding.

    Metrics computed:
      Batch correction : iLISI, graph_connectivity
      Bio conservation : cLISI, ASW_label (when cell_type_key is provided)

    kBET is not computed by default — it can be very slow on large datasets.

    Parameters
    ----------
    mdata : MuData
        MuData with the embedding in ``obsm[embed_key]`` and batch/cell-type
        columns accessible from ``mdata.obs`` or ``mdata["rna"].obs``.
    embed_key : str
        Key in ``mdata.obsm`` holding the embedding to benchmark.
    batch_key : str
        obs column for batch labels.
    cell_type_key : str or None
        obs column for cell type labels.  Required for cLISI and ASW_label.
        Pass ``None`` to skip bio-conservation metrics.
    label : str
        Short identifier added to all metric names (e.g. "mofa", "totalvi").

    Returns
    -------
    dict
        Keys: ``{label}_ilisi``, ``{label}_graph_conn``,
              ``{label}_clisi`` (if cell_type_key),
              ``{label}_asw_label`` (if cell_type_key),
              ``{label}_scib_error`` (if scib raises).
    """
    try:
        import scib
    except ImportError:
        return {f"{label}_scib_error": "scib not installed (pip install scib)"}

    if embed_key not in mdata.obsm:
        return {f"{label}_scib_error": f"{embed_key} not found in mdata.obsm"}

    # Build a temporary AnnData that scib expects:
    #   .X or .obsm[embed_key]  — embedding
    #   .obs[batch_key]         — batch labels
    #   .obs[cell_type_key]     — cell type labels (optional)
    import anndata as ad

    obs = mdata.obs.copy() if batch_key in mdata.obs.columns \
        else mdata["rna"].obs.copy()

    adata_tmp = ad.AnnData(
        X=mdata.obsm[embed_key].copy(),
        obs=obs.loc[mdata.obs_names] if set(mdata.obs_names).issubset(obs.index)
            else obs,
    )
    adata_tmp.obsm[embed_key] = mdata.obsm[embed_key].copy()

    # Ensure batch_key is present
    if batch_key not in adata_tmp.obs.columns:
        rna_obs = mdata["rna"].obs
        if batch_key in rna_obs.columns:
            adata_tmp.obs[batch_key] = rna_obs.loc[
                adata_tmp.obs_names, batch_key
            ].values
        else:
            return {f"{label}_scib_error": f"batch_key='{batch_key}' not found"}

    # scib internally calls .cat.categories / .value_counts() on label columns,
    # which requires CategoricalDtype — plain str (object) raises
    # "Can only use .cat accessor with a 'category' dtype".
    # Cast every label column to Categorical immediately after building adata_tmp.
    import pandas as pd
    adata_tmp.obs[batch_key] = pd.Categorical(
        adata_tmp.obs[batch_key].astype(str)
    )

    # Ensure cell_type_key is present if requested
    if cell_type_key is not None and cell_type_key not in adata_tmp.obs.columns:
        for src_obs in [mdata.obs, mdata["rna"].obs, mdata["adt"].obs]:
            if cell_type_key in src_obs.columns:
                adata_tmp.obs[cell_type_key] = src_obs.loc[
                    adata_tmp.obs_names, cell_type_key
                ].values
                break
        else:
            cell_type_key = None  # silently skip bio metrics

    # Build neighbor graph on the embedding (required by scib)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(adata_tmp, use_rep=embed_key, random_state=0)

    # Wrap every scib call in catch_warnings — scib internals emit FutureWarnings
    # (e.g. pd.value_counts deprecated in graph_connectivity.py) that escape the
    # runner's global filterwarnings because they originate inside scib's own module.
    results: dict = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        try:
            # ── Batch correction ──
            results[f"{label}_ilisi"] = float(
                np.mean(scib.metrics.ilisi_graph(
                    adata_tmp, batch_key=batch_key, type_="embed",
                    use_rep=embed_key,
                ))
            )
        except Exception as exc:
            results[f"{label}_ilisi_error"] = str(exc)

        try:
            results[f"{label}_graph_conn"] = float(
                scib.metrics.graph_connectivity(adata_tmp, label_key=batch_key)
            )
        except Exception as exc:
            results[f"{label}_graph_conn_error"] = str(exc)

        if cell_type_key is not None:
            # Also cast cell_type_key to Categorical for clisi_graph / silhouette
            adata_tmp.obs[cell_type_key] = pd.Categorical(
                adata_tmp.obs[cell_type_key].astype(str)
            )

            try:
                results[f"{label}_clisi"] = float(
                    np.mean(scib.metrics.clisi_graph(
                        adata_tmp, label_key=cell_type_key, type_="embed",
                        use_rep=embed_key,
                    ))
                )
            except Exception as exc:
                results[f"{label}_clisi_error"] = str(exc)

            try:
                results[f"{label}_asw_label"] = float(
                    scib.metrics.silhouette(
                        adata_tmp, label_key=cell_type_key, embed=embed_key,
                    )
                )
            except Exception as exc:
                results[f"{label}_asw_label_error"] = str(exc)

    return results


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
    compute_scib: bool = False,
    cell_type_key: Optional[str] = None,
) -> tuple[MuData, dict]:
    """
    Run MOFA+ integration on paired RNA + ADT data.

    Parameters
    ----------
    mdata : MuData
        MuData with ``mdata["rna"]`` and ``mdata["adt"]`` modalities.
        - ``mdata["rna"].X``     — log1p-normalized RNA values
        - ``mdata["adt"].X``     — CLR-normalized ADT values
        - ``obs[batch_key]``     — batch column (present on both modalities)
    batch_key : str
        Column in ``.obs`` used as the group/batch variable for MOFA+.
    n_factors : int
        Number of MOFA+ latent factors.  Default: 15.
    use_layer : str or None
        Layer to use instead of ``.X``.  Default: None.
    random_state : int
        Random seed.  Default: 0.
    inplace : bool
        Modify ``mdata`` in place if True.
    compute_scib : bool
        Compute scib benchmark metrics after integration.  Default: False.
    cell_type_key : str or None
        obs column for cell types — required for bio-conservation scib metrics.

    Returns
    -------
    mdata_out : MuData
        With ``obsm["X_mofa"]``, ``obsm["X_umap_mofa"]``,
        ``uns["omicsage_mofa"]``.
    metrics : dict
        Includes ``metrics["scib"]`` dict when ``compute_scib=True``.
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(mdata_out, use_rep="X_mofa", random_state=random_state)
        sc.tl.umap(mdata_out, random_state=random_state)

    mdata_out.obsm["X_umap_mofa"] = mdata_out.obsm.pop("X_umap")

    if "LFs" in mdata_out.varm:
        print(f"[run_mofa] varm['LFs'] written: shape {mdata_out.varm['LFs'].shape}")
    else:
        print("[run_mofa] WARNING: mdata.varm['LFs'] missing — trying fallback")
        rna_lfs = mdata_out["rna"].varm.get("LFs")
        adt_lfs = mdata_out["adt"].varm.get("LFs")
        if rna_lfs is not None and adt_lfs is not None:
            mdata_out.varm["LFs"] = np.concatenate(
                [np.asarray(rna_lfs), np.asarray(adt_lfs)], axis=0
            )
            print(f"[run_mofa] Built mdata.varm['LFs'] from modality varms: "
                  f"shape {mdata_out.varm['LFs'].shape}")
        else:
            print("[run_mofa] WARNING: LFs not found anywhere — "
                  "weights heatmap will show placeholder in report")

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
        "method": "mofa",
    }

    if compute_scib:
        print("[run_mofa] Computing scib metrics ...", flush=True)
        scib_scores = _compute_scib_metrics(
            mdata_out, embed_key="X_mofa",
            batch_key=batch_key, cell_type_key=cell_type_key,
            label="mofa",
        )
        metrics["scib"] = scib_scores
        print(f"[run_mofa] scib: {scib_scores}", flush=True)

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
    compute_scib: bool = False,
    cell_type_key: Optional[str] = None,
) -> tuple[MuData, dict]:
    """
    Run totalVI integration on paired RNA + ADT data.

    Parameters
    ----------
    mdata : MuData
        MuData with ``mdata["rna"]`` and ``mdata["adt"]`` modalities.
        - ``mdata["rna"].layers["counts"]`` — raw RNA integer counts (required)
        - ``mdata["adt"].layers["counts"]`` — raw ADT integer counts (required)
        - ``mdata["rna"].obs[batch_key]``   — batch column on RNA modality
    batch_key : str
        Column in ``mdata["rna"].obs`` used as batch variable.
    max_epochs : int
        Maximum training epochs.  Default: 400.
    random_state : int
        Random seed for scvi-tools.  Default: 0.
    inplace : bool
        Modify ``mdata`` in place if True.
    compute_scib : bool
        Compute scib benchmark metrics after integration.  Default: False.
    cell_type_key : str or None
        obs column for cell types — required for bio-conservation scib metrics.

    Returns
    -------
    mdata_out : MuData
        With ``obsm["X_totalVI"]``, ``obsm["X_umap_totalVI"]``,
        ``uns["omicsage_totalVI"]``.
    metrics : dict
        Includes ``metrics["scib"]`` dict when ``compute_scib=True``.
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
        "method": "totalvi",
    }

    if compute_scib:
        print("[run_totalvi] Computing scib metrics ...", flush=True)
        scib_scores = _compute_scib_metrics(
            mdata_out, embed_key="X_totalVI",
            batch_key=batch_key, cell_type_key=cell_type_key,
            label="totalvi",
        )
        metrics["scib"] = scib_scores
        print(f"[run_totalvi] scib: {scib_scores}", flush=True)

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
# C. Both (MOFA+ then totalVI)
# ---------------------------------------------------------------------------

def run_both(
    mdata: MuData,
    batch_key: str,
    n_factors: int = 15,
    max_epochs: int = 400,
    random_state: int = 0,
    inplace: bool = False,
    compute_scib: bool = False,
    cell_type_key: Optional[str] = None,
) -> tuple[MuData, dict]:
    """
    Run MOFA+ **and** totalVI sequentially on the same MuData.

    Both embeddings and their UMAPs are written.  When ``compute_scib=True``,
    scib metrics are computed for both methods and a side-by-side comparison
    is included in ``metrics["scib_comparison"]``.

    Parameters
    ----------
    mdata : MuData
        Same prerequisites as both ``run_mofa`` and ``run_totalvi``:
        needs both log1p-normalised X AND raw counts layers.
    batch_key : str
        obs column for batch.
    n_factors : int
        MOFA+ latent factors.  Default: 15.
    max_epochs : int
        totalVI training epochs.  Default: 400.
    random_state : int
        Random seed for both methods.  Default: 0.
    inplace : bool
        Modify ``mdata`` in place if True.
    compute_scib : bool
        Compute scib metrics for both embeddings.  Default: False.
    cell_type_key : str or None
        obs column for cell types (scib bio-conservation metrics).

    Returns
    -------
    mdata_out : MuData
        Contains all four embeddings:
          ``obsm["X_mofa"]``         ``obsm["X_umap_mofa"]``
          ``obsm["X_totalVI"]``      ``obsm["X_umap_totalVI"]``
    metrics : dict
        ``metrics["mofa"]``   — MOFA+ metrics dict
        ``metrics["totalvi"]`` — totalVI metrics dict
        ``metrics["method"]``  == "both"
        ``metrics["scib_comparison"]`` — flat dict merging both scib dicts
                                          (only when compute_scib=True)
    """
    print("[run_both] Running MOFA+ ...", flush=True)
    mdata_out, mofa_metrics = run_mofa(
        mdata,
        batch_key=batch_key,
        n_factors=n_factors,
        random_state=random_state,
        inplace=inplace,
        compute_scib=compute_scib,
        cell_type_key=cell_type_key,
    )

    print("[run_both] Running totalVI ...", flush=True)
    mdata_out, totalvi_metrics = run_totalvi(
        mdata_out,
        batch_key=batch_key,
        max_epochs=max_epochs,
        random_state=random_state,
        inplace=True,   # already operating on the copy from run_mofa
        compute_scib=compute_scib,
        cell_type_key=cell_type_key,
    )

    combined_metrics: dict = {
        "method": "both",
        "batch_key": batch_key,
        "n_cells": int(mdata_out.n_obs),
        "mofa": mofa_metrics,
        "totalvi": totalvi_metrics,
    }

    if compute_scib:
        mofa_scib   = mofa_metrics.get("scib", {})
        totalvi_scib = totalvi_metrics.get("scib", {})
        combined_metrics["scib_comparison"] = {**mofa_scib, **totalvi_scib}
        print("[run_both] scib comparison:", combined_metrics["scib_comparison"],
              flush=True)

    return mdata_out, combined_metrics


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
