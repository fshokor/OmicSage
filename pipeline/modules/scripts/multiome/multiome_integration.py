"""
multiome_integration.py — Multi-modal RNA + ATAC integration for OmicSage.

Implements two fully-Python integration methods for paired RNA + ATAC multiome data:

  A. MOFA+  (Multi-Omics Factor Analysis)
     Linear factor model. Fast, interpretable, handles batch via groups_label.
     Reference: Argelaguet et al., Genome Biology 2020.
     Python: mu.tl.mofa(mdata, groups_label=batch_key)
     Input:  mdata["rna"].X (log1p-normalised), mdata["atac"].X (TF-IDF)

  B. MultiVI (Multi-modal Variational Inference)
     Deep generative model for paired RNA + ATAC.
     Models RNA with NB distribution, ATAC with Bernoulli (binarised peaks).
     Requires raw counts for RNA (NB) and raw counts for ATAC (binarised).
     Reference: Ashuach et al., Nature Methods 2023.
     Python: scvi.model.MULTIVI

Embedding keys written to mdata (OmicSage convention):
  mdata.obsm["X_mofa"]          — MOFA+ latent factors
  mdata.obsm["X_umap_mofa"]     — UMAP from MOFA+ embedding
  mdata.obsm["X_multivi"]       — MultiVI latent representation
  mdata.obsm["X_umap_multivi"]  — UMAP from MultiVI embedding

API
---
run_mofa(mdata, batch_key, n_factors=15, use_layer=None, random_state=0,
         inplace=False) -> (MuData, dict)

run_multivi(mdata, batch_key, n_latent=20, max_epochs=500, random_state=0,
            inplace=False) -> (MuData, dict)

References
----------
sc-best-practices paired integration (MultiVI section):
  https://www.sc-best-practices.org/multimodal_integration/paired_integration.html#multiome-data
MultiVI paper:
  Ashuach et al., Nature Methods 2023
  https://www.nature.com/articles/s41592-023-02011-6
scvi-tools MultiVI tutorial:
  https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/MultiVI_tutorial.html
MOFA+ paper:
  Argelaguet et al., Genome Biology 2020
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
    _MULTIVI = _scvi_model.MULTIVI
    _SCVI_AVAILABLE = True
except (ModuleNotFoundError, ImportError, Exception):
    _SCVI_AVAILABLE = False
    _MULTIVI = None


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
    Run MOFA+ integration on paired RNA + ATAC multiome data.

    MOFA+ learns a set of latent factors that explain variance across both
    modalities jointly.  The batch_key is passed as groups_label so MOFA+
    treats each batch as a separate group, enabling batch-aware factorisation.

    Parameters
    ----------
    mdata : MuData
        MuData with ``mdata["rna"]`` and ``mdata["atac"]`` modalities.
        - ``mdata["rna"].X``     — log1p-normalised RNA values
        - ``mdata["atac"].X``    — TF-IDF normalised ATAC values
        - ``obs[batch_key]``     — batch column (present on RNA modality at minimum)
    batch_key : str
        Column in ``mdata["rna"].obs`` used as the group/batch variable for MOFA+.
    n_factors : int
        Number of MOFA+ latent factors.  Default: 15.
    use_layer : str or None
        Layer to use instead of ``.X``.  Default: None (use .X).
    random_state : int
        Random seed.  Default: 0.
    inplace : bool
        Modify ``mdata`` in place if True; default makes a copy.

    Returns
    -------
    mdata_out : MuData
        With ``obsm["X_mofa"]``, ``obsm["X_umap_mofa"]``,
        ``uns["omicsage_mofa"]``.
    metrics : dict
        n_cells, n_factors, batch_key, n_batches, batch_values,
        mofa_key, umap_key, umap_computed, random_state, method.

    Raises
    ------
    ImportError
        If muon or mofapy2 are not installed.
    KeyError
        If "rna" or "atac" modalities are missing, or batch_key not in obs.
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

    # Push batch key to top-level mdata.obs — muon namespaces modality obs
    # columns as "rna:batch" / "atac:batch" which MOFA+ cannot see directly.
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

    # Rename X_umap → X_umap_mofa to avoid collision with RNA/ATAC UMAPs
    mdata_out.obsm["X_umap_mofa"] = mdata_out.obsm.pop("X_umap")

    batch_values = sorted(
        mdata_out["rna"].obs[batch_key].astype(str).unique().tolist()
    )
    n_factors_actual = mdata_out.obsm["X_mofa"].shape[1]

    metrics: dict = {
        "n_cells":       int(mdata_out.n_obs),
        "n_factors":     n_factors_actual,
        "batch_key":     batch_key,
        "n_batches":     len(batch_values),
        "batch_values":  batch_values,
        "mofa_key":      "X_mofa",
        "umap_key":      "X_umap_mofa",
        "umap_computed": True,
        "random_state":  random_state,
        "method":        "mofa",
    }

    mdata_out.uns["omicsage_mofa"] = {
        "module":    "multiome_integration.run_mofa",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "batch_key":    batch_key,
            "n_factors":    n_factors,
            "use_layer":    use_layer,
            "random_state": random_state,
        },
        "outputs": {
            "mofa_key":        "X_mofa",
            "umap_key":        "X_umap_mofa",
            "n_factors_actual": n_factors_actual,
        },
    }

    return mdata_out, metrics


# ---------------------------------------------------------------------------
# B. MultiVI
# ---------------------------------------------------------------------------

def run_multivi(
    mdata: MuData,
    batch_key: str,
    n_latent: int = 20,
    max_epochs: int = 500,
    random_state: int = 0,
    inplace: bool = False,
) -> tuple[MuData, dict]:
    """
    Run MultiVI integration on paired RNA + ATAC multiome data.

    MultiVI learns a joint latent space by modelling RNA with a
    negative-binomial distribution and ATAC accessibility with a Bernoulli
    distribution (after binarisation).  It is the recommended deep generative
    model for paired multiome data per sc-best-practices.

    Implementation approach
    -----------------------
    MultiVI is set up via ``scvi.model.MULTIVI.setup_mudata`` which accepts
    a MuData directly, designating which modality provides RNA counts and
    which provides ATAC counts.  This avoids the need to concatenate RNA
    genes + ATAC peaks into a single AnnData.

    ATAC peaks are binarised internally by MultiVI (any count > 0 → 1).
    Raw integer counts must be present in:
      mdata["rna"].layers["counts"]   — raw RNA counts (not log-normalised)
      mdata["atac"].layers["counts"]  — raw ATAC peak counts

    Parameters
    ----------
    mdata : MuData
        MuData with ``mdata["rna"]`` and ``mdata["atac"]`` modalities.
        - ``mdata["rna"].layers["counts"]``  — raw RNA integer counts (required)
        - ``mdata["atac"].layers["counts"]`` — raw ATAC integer counts (required)
        - ``mdata["rna"].obs[batch_key]``    — batch column on RNA modality
    batch_key : str
        Column in ``mdata["rna"].obs`` used as batch variable.
    n_latent : int
        Dimensionality of the latent space.  Default: 20.
    max_epochs : int
        Maximum training epochs.  Default: 500.
    random_state : int
        Random seed for scvi-tools.  Default: 0.
    inplace : bool
        Modify ``mdata`` in place if True; default makes a copy.

    Returns
    -------
    mdata_out : MuData
        With ``obsm["X_multivi"]``, ``obsm["X_umap_multivi"]``,
        ``uns["omicsage_multivi"]``.
    metrics : dict
        n_cells, n_latent, batch_key, n_batches, batch_values,
        max_epochs, multivi_key, umap_key, umap_computed,
        random_state, method.

    Raises
    ------
    ImportError
        If scvi-tools is not installed.
    KeyError
        If "rna" or "atac" modalities are missing, or required layers/keys absent.
    """
    _validate_multivi_inputs(mdata, batch_key)

    try:
        import scvi
        from scvi.model import MULTIVI
    except (ModuleNotFoundError, ImportError) as exc:
        raise ImportError(
            "scvi-tools is required for MultiVI integration. "
            "Install with: pip install scvi-tools"
        ) from exc

    mdata_out = mdata if inplace else mdata.copy()

    rna  = mdata_out["rna"]
    atac = mdata_out["atac"]

    if hasattr(scvi, "settings"):
        scvi.settings.seed = random_state

    # Push batch key to top-level mdata.obs (same pattern as MOFA+)
    mdata_out.obs[batch_key] = rna.obs[batch_key].values

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # setup_mudata registers the MuData with MultiVI:
        #   rna_layer="counts"       → raw RNA counts for NB model
        #   atac_layer="counts"      → raw ATAC counts, binarised internally
        #   batch_key from rna obs   → batch correction
        MULTIVI.setup_mudata(
            mdata_out,
            rna_layer="counts",
            atac_layer="counts",
            batch_key=batch_key,
            modalities={
                "rna_layer":   "rna",
                "batch_key":   "rna",
                "atac_layer":  "atac",
            },
        )

        model = MULTIVI(
            mdata_out,
            n_latent=n_latent,
        )
        model.train(max_epochs=max_epochs)

    latent = model.get_latent_representation()  # (n_cells, n_latent)
    mdata_out.obsm["X_multivi"] = latent

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(mdata_out, use_rep="X_multivi", random_state=random_state)
        sc.tl.umap(mdata_out, random_state=random_state)

    mdata_out.obsm["X_umap_multivi"] = mdata_out.obsm.pop("X_umap")

    batch_values = sorted(rna.obs[batch_key].astype(str).unique().tolist())
    n_latent_actual = latent.shape[1]

    metrics: dict = {
        "n_cells":       int(mdata_out.n_obs),
        "n_latent":      n_latent_actual,
        "batch_key":     batch_key,
        "n_batches":     len(batch_values),
        "batch_values":  batch_values,
        "max_epochs":    max_epochs,
        "multivi_key":   "X_multivi",
        "umap_key":      "X_umap_multivi",
        "umap_computed": True,
        "random_state":  random_state,
        "method":        "multivi",
    }

    mdata_out.uns["omicsage_multivi"] = {
        "module":    "multiome_integration.run_multivi",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "batch_key":    batch_key,
            "n_latent":     n_latent,
            "max_epochs":   max_epochs,
            "random_state": random_state,
        },
        "outputs": {
            "multivi_key":      "X_multivi",
            "umap_key":         "X_umap_multivi",
            "n_latent_actual":  n_latent_actual,
        },
    }

    return mdata_out, metrics


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def _validate_mofa_inputs(mdata: MuData, batch_key: str) -> None:
    for mod in ("rna", "atac"):
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


def _validate_multivi_inputs(mdata: MuData, batch_key: str) -> None:
    for mod in ("rna", "atac"):
        if mod not in mdata.mod:
            raise KeyError(
                f"mdata must contain a '{mod}' modality. "
                f"Found: {list(mdata.mod.keys())}"
            )
    if "counts" not in mdata["rna"].layers:
        raise KeyError(
            "mdata['rna'].layers['counts'] not found. "
            "MultiVI requires raw RNA counts. "
            "Ensure ingest.py preserved raw counts in layers['counts']."
        )
    if "counts" not in mdata["atac"].layers:
        raise KeyError(
            "mdata['atac'].layers['counts'] not found. "
            "MultiVI requires raw ATAC peak counts. "
            "Ensure atac_qc.py preserved raw counts in layers['counts']."
        )
    if batch_key not in mdata["rna"].obs.columns:
        raise KeyError(
            f"batch_key='{batch_key}' not found in mdata['rna'].obs. "
            f"Available: {list(mdata['rna'].obs.columns)}"
        )
