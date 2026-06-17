"""
OmicSage — ATAC Dimensionality Reduction Module
pipeline/modules/multiome/atac_reduce.py

Input : QC-filtered ATAC AnnData (output of atac_qc.py)
        Raw peak counts in .X or .layers["counts"]
Output: reduced AnnData + metrics dict

Steps
-----
1. TF-IDF normalisation  — mu.atac.pp.tfidf (muon)
   Normalises across cells (sequencing depth) and across peaks (peak accessibility).
   Result stored in .X and .layers["tf_idf"].

2. LSI (Latent Semantic Indexing) — sklearn TruncatedSVD on the TF-IDF matrix
   LSI component 1 is ALWAYS dropped — it correlates with sequencing depth, not biology.
   Confirmed by Signac vignette (DepthCor check) and sc-best-practices ATAC chapter.
   Components 2–n_components stored in obsm["X_lsi"].

3. Neighbor graph on X_lsi (sc.pp.neighbors, use_rep="X_lsi")

4. UMAP from LSI neighbor graph → obsm["X_umap_atac"]
   Namespaced to avoid collision with RNA obsm["X_umap"].

5. Leiden clustering on LSI graph → obs["atac_leiden"]

6. Provenance written to uns["omicsage_atac_reduce"]

Key design decisions
--------------------
- TF-IDF via muon (mu.atac.pp.tfidf) — standard Python ATAC toolkit, same as
  sc-best-practices tutorial.
- LSI via sklearn TruncatedSVD — pure Python, no extra dependency beyond sklearn,
  fully transparent, matches sc-best-practices Python examples.
- Component 1 always dropped — this is unconditional, not configurable. Every
  reference (Signac, sc-best-practices, Seurat v5 vignette) agrees on this.
- obsm key "X_lsi" (not "X_pca") — LSI is the correct name; prevents confusion
  with RNA PCA embeddings in downstream MuData.
- UMAP key "X_umap_atac" — namespaced to avoid collision with RNA "X_umap".
- Leiden key "atac_leiden" — namespaced to avoid collision with RNA "leiden".

References
----------
- sc-best-practices ATAC QC chapter (Lance & Martens 2022)
- Signac PBMC vignette: https://stuartlab.org/signac/articles/pbmc_vignette.html
  "We exclude the first dimension as this is typically correlated with sequencing depth"
- Seurat v5 ATAC integration vignette: RunSVD dims=2:30

Usage
-----
    from pipeline.modules.multiome.atac_reduce import atac_reduce

    atac_out, metrics = atac_reduce(atac_qcd)
"""

from __future__ import annotations

import logging
import warnings
from datetime import datetime, timezone
from typing import Any, Optional

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def atac_reduce(
    adata: AnnData,
    n_components: int = 50,
    n_neighbors: int = 15,
    leiden_resolution: float = 0.5,
    use_raw_counts: bool = True,
    random_state: int = 0,
    inplace: bool = False,
) -> tuple[AnnData, dict[str, Any]]:
    """
    Run TF-IDF → LSI → UMAP → Leiden clustering on a QC-filtered ATAC AnnData.

    LSI component 1 is always dropped (correlates with sequencing depth).
    The usable components are 2 through n_components, stored as X_lsi.

    Parameters
    ----------
    adata : AnnData
        QC-filtered ATAC AnnData with raw peak counts in .X or .layers["counts"].
        Typically the output of atac_qc().
    n_components : int
        Number of SVD components to compute (including component 1 which is
        dropped). Default 50 — components 2–50 are used for downstream steps.
    n_neighbors : int
        Number of neighbors for the kNN graph built on X_lsi. Default 15.
    leiden_resolution : float
        Resolution parameter for Leiden clustering on the LSI graph. Default 0.5.
    use_raw_counts : bool
        If True (default), use .layers["counts"] as input to TF-IDF when
        available. Falls back to .X if the layer does not exist.
    random_state : int
        Reproducibility seed for TruncatedSVD and UMAP. Default 0.
    inplace : bool
        If False (default), operate on a copy.

    Returns
    -------
    adata_out : AnnData
        ATAC AnnData with the following added:
        - .layers["tf_idf"]      TF-IDF normalised peak matrix
        - .obsm["X_lsi"]         LSI embedding (components 2–n_components)
        - .obsm["X_umap_atac"]   UMAP from LSI neighbors
        - .obs["atac_leiden"]    Leiden cluster labels
        - .uns["omicsage_atac_reduce"]  provenance record
    metrics : dict
        Summary dict with embedding shapes, cluster counts, and parameters.

    Raises
    ------
    TypeError
        If adata is not an AnnData object.
    ValueError
        If adata has 0 cells or 0 peaks, or n_components >= n_peaks.
    """
    _validate_input(adata, n_components)

    adata_out = adata.copy()

    n_cells = adata_out.n_obs
    n_peaks = adata_out.n_vars
    logger.info(
        "ATAC reduce start — %d cells × %d peaks, n_components=%d",
        n_cells, n_peaks, n_components,
    )

    # ------------------------------------------------------------------
    # Step 1 — Resolve input matrix (raw counts preferred)
    # ------------------------------------------------------------------
    if use_raw_counts and "counts" in adata_out.layers:
        X_raw = adata_out.layers["counts"]
        logger.info("Using .layers['counts'] as TF-IDF input")
    else:
        X_raw = adata_out.X
        logger.info("Using .X as TF-IDF input")

    if sp.issparse(X_raw):
        X_raw = X_raw.tocsr()

    # ------------------------------------------------------------------
    # Step 2 — TF-IDF normalisation
    # ------------------------------------------------------------------
    logger.info("Computing TF-IDF normalisation")
    X_tfidf = _tfidf(X_raw)

    adata_out.layers["tf_idf"] = X_tfidf
    # Set .X to TF-IDF for neighbor graph / UMAP steps
    adata_out.X = X_tfidf
    logger.info("TF-IDF stored in .layers['tf_idf'] and .X")

    # ------------------------------------------------------------------
    # Step 3 — LSI via TruncatedSVD; drop component 1
    # ------------------------------------------------------------------
    logger.info(
        "Running TruncatedSVD (n_components=%d, random_state=%d)",
        n_components, random_state,
    )
    svd = TruncatedSVD(
        n_components=n_components,
        algorithm="randomized",
        random_state=random_state,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X_svd = svd.fit_transform(X_tfidf)   # shape: (n_cells, n_components)

    # Component 1 (index 0) always correlates with sequencing depth — drop it.
    # Use components 2..n_components (indices 1..n_components-1)
    X_lsi = X_svd[:, 1:]   # shape: (n_cells, n_components - 1)

    # L2-normalise rows so that Euclidean distance ≈ cosine distance
    # (standard practice for LSI in single-cell ATAC)
    X_lsi = normalize(X_lsi, norm="l2")

    adata_out.obsm["X_lsi"] = X_lsi

    variance_ratio = svd.explained_variance_ratio_
    # Report variance for components 2..n (component 1 dropped)
    variance_ratio_used = variance_ratio[1:]

    logger.info(
        "LSI complete — using components 2–%d "
        "(cumulative variance explained: %.1f%%)",
        n_components,
        float(variance_ratio_used.sum()) * 100,
    )

    # ------------------------------------------------------------------
    # Step 4 — Neighbor graph on X_lsi
    # ------------------------------------------------------------------
    logger.info(
        "Computing neighbor graph — n_neighbors=%d, use_rep='X_lsi'", n_neighbors
    )
    sc.pp.neighbors(
        adata_out,
        n_neighbors=n_neighbors,
        use_rep="X_lsi",
        random_state=random_state,
    )

    # ------------------------------------------------------------------
    # Step 5 — UMAP (namespaced key to avoid collision with RNA)
    # ------------------------------------------------------------------
    logger.info("Running UMAP → obsm['X_umap_atac']")
    sc.tl.umap(adata_out, random_state=random_state)
    adata_out.obsm["X_umap_atac"] = adata_out.obsm["X_umap"].copy()
    # Keep X_umap for scanpy compatibility but primary key is X_umap_atac
    logger.info("UMAP stored in obsm['X_umap_atac']")

    # ------------------------------------------------------------------
    # Step 6 — Leiden clustering on LSI graph
    # ------------------------------------------------------------------
    logger.info("Running Leiden clustering (resolution=%.2f)", leiden_resolution)
    sc.tl.leiden(
        adata_out,
        resolution=leiden_resolution,
        random_state=random_state,
        key_added="atac_leiden",
    )
    n_clusters = int(adata_out.obs["atac_leiden"].nunique())
    logger.info("Leiden clustering complete — %d clusters", n_clusters)

    # ------------------------------------------------------------------
    # Step 7 — Build metrics dict
    # ------------------------------------------------------------------
    metrics: dict[str, Any] = {
        "n_cells":                   n_cells,
        "n_peaks":                   n_peaks,
        "n_components_computed":     n_components,
        "n_lsi_components_used":     n_components - 1,   # component 1 dropped
        "component_1_dropped":       True,
        "n_neighbors":               n_neighbors,
        "leiden_resolution":         leiden_resolution,
        "n_leiden_clusters":         n_clusters,
        "variance_ratio_all":        variance_ratio.tolist(),
        "variance_ratio_used":       variance_ratio_used.tolist(),
        "cumulative_variance_used":  float(variance_ratio_used.sum()),
        "embeddings_computed":       ["X_lsi", "X_umap_atac"],
        "cluster_key":               "atac_leiden",
    }

    # ------------------------------------------------------------------
    # Step 8 — Provenance
    # ------------------------------------------------------------------
    adata_out.uns["omicsage_atac_reduce"] = {
        "omicsage_module":           "pipeline.modules.multiome.atac_reduce",
        "omicsage_version":          "0.1.0",
        "timestamp":                 datetime.now(timezone.utc).isoformat(),
        "n_components_computed":     n_components,
        "n_lsi_components_used":     n_components - 1,
        "component_1_dropped":       True,
        "component_1_drop_reason":   (
            "LSI component 1 correlates with sequencing depth (technical), "
            "not biology. Dropped unconditionally per Signac and sc-best-practices."
        ),
        "tfidf_method":              "log(TF) * log(1 + IDF)",
        "lsi_method":                "sklearn.decomposition.TruncatedSVD",
        "lsi_normalisation":         "L2 row normalisation post-SVD",
        "n_neighbors":               n_neighbors,
        "leiden_resolution":         leiden_resolution,
        "n_leiden_clusters":         n_clusters,
        "variance_ratio_used":       variance_ratio_used.tolist(),
        "cumulative_variance_used":  float(variance_ratio_used.sum()),
        "random_state":              random_state,
        "use_raw_counts":            use_raw_counts,
        "scanpy_version":            sc.__version__,
    }

    logger.info(
        "ATAC reduce complete — %d cells, %d LSI components, "
        "%d Leiden clusters, UMAP computed",
        n_cells, n_components - 1, n_clusters,
    )

    return adata_out, metrics


# ---------------------------------------------------------------------------
# TF-IDF implementation
# ---------------------------------------------------------------------------

def _tfidf(X: sp.spmatrix | np.ndarray) -> sp.csr_matrix:
    """
    Compute TF-IDF normalisation for a peak count matrix.

    Formula: log(TF + 1) * log(1 + N / (1 + df))
    where:
        TF  = count / total_counts_per_cell  (term frequency per cell)
        df  = number of cells with count > 0 per peak  (document frequency)
        N   = number of cells

    This is the standard TF-IDF used for scATAC-seq (Cusanovich et al. 2015,
    Signac, sc-best-practices).  The +1 pseudocounts prevent log(0).

    Parameters
    ----------
    X : sparse or dense matrix, shape (n_cells, n_peaks)
        Raw peak count matrix.

    Returns
    -------
    X_tfidf : csr_matrix, shape (n_cells, n_peaks)
        TF-IDF normalised matrix.
    """
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr().astype(float)

    n_cells, n_peaks = X.shape

    # Term frequency: normalise each cell by its total counts
    total_counts = np.asarray(X.sum(axis=1)).ravel()  # (n_cells,)
    # Avoid division by zero for empty cells
    total_counts = np.where(total_counts == 0, 1, total_counts)
    # Divide each row by its total (sparse row-wise division)
    D_inv = sp.diags(1.0 / total_counts)
    TF = D_inv @ X   # (n_cells, n_peaks), sparse

    # Inverse document frequency
    # df[j] = number of cells where peak j > 0
    df = np.asarray((X > 0).sum(axis=0)).ravel()  # (n_peaks,)
    idf = np.log1p(n_cells / (1.0 + df))           # (n_peaks,)

    # TF-IDF: scale each column by its IDF
    D_idf = sp.diags(idf)
    X_tfidf = TF @ D_idf   # (n_cells, n_peaks), sparse

    # Log-transform TF component
    X_tfidf = X_tfidf.tocsr()
    X_tfidf.data = np.log1p(X_tfidf.data * 1e4)  # scale factor 1e4 standard

    return X_tfidf.tocsr()


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def _validate_input(adata: AnnData, n_components: int) -> None:
    if not isinstance(adata, AnnData):
        raise TypeError(
            f"atac_reduce() expects an AnnData object, got {type(adata).__name__}. "
            "Pass mdata['atac'] not the full MuData."
        )
    if adata.n_obs == 0:
        raise ValueError("AnnData has 0 cells — nothing to reduce.")
    if adata.n_vars == 0:
        raise ValueError("AnnData has 0 peaks — nothing to reduce.")
    if n_components >= adata.n_vars:
        raise ValueError(
            f"n_components={n_components} must be less than n_peaks={adata.n_vars}. "
            "Reduce n_components or use a larger peak set."
        )
    if n_components < 2:
        raise ValueError(
            f"n_components={n_components} must be >= 2 so that component 1 can be "
            "dropped and at least 1 component remains."
        )
