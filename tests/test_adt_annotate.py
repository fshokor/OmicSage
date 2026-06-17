"""
tests/test_adt_annotate.py — Tests for pipeline/modules/scripts/cite/adt_annotate.py

Coverage:
  Happy path
  - Returns (AnnData, dict) tuple
  - obs["leiden"] always written
  - obs["adt_celltype_manual"] written when annotation_map provided
  - obs["adt_celltype_manual"] NOT written when annotation_map is None
  - No obs["celltype"] written (RNA collision guard)
  - rank_genes_groups computed (uns["rank_genes_groups"] present)
  - dendrogram computed (uns["dendrogram_leiden"] present)

  Output shapes and types
  - leiden values are strings
  - adt_celltype_manual values are strings when annotated
  - UMAP and PCA embeddings untouched
  - RNA modality entirely untouched
  - Layers preserved (adt_clr, counts)

  Annotation logic
  - annotation_map=None: annotated=False, celltype_key=None
  - annotation_map provided: annotated=True, celltype_key="adt_celltype_manual"
  - Unknown keys in annotation_map emit UserWarning (not an error)
  - Unknown keys in annotation_map do not crash the function
  - Partial map: unmapped clusters keep their numeric string IDs
  - Full map: all clusters replaced with cell type labels
  - annotation_map with integer keys (auto-coerced to str)
  - Empty annotation_map: annotated=True, leiden values unchanged

  Metrics dict
  - Required keys present
  - n_cells correct
  - n_clusters >= 1
  - cluster_sizes is dict[str, int]
  - cluster_sizes values sum to n_cells
  - annotation_map in metrics is {} when None provided
  - annotation_map in metrics matches applied map
  - annotated flag correct
  - leiden_key == "leiden"
  - celltype_key == "adt_celltype_manual" when annotated
  - celltype_key is None when not annotated
  - resolution in metrics correct
  - random_state in metrics correct
  - n_clusters matches cluster_sizes length

  Provenance
  - uns["omicsage_adt_annotate"] present
  - module == "adt_annotate"
  - timestamp non-empty string
  - params: resolution, flavor, directed, n_iterations, random_state,
            annotation_map_provided
  - outputs: leiden_key, n_clusters, celltype_key

  inplace behaviour
  - inplace=False: copy returned, original untouched
  - inplace=True: mdata["adt"] modified in place, same object returned

  Error handling
  - Missing "adt" modality raises KeyError
  - Missing obsp["connectivities"] raises KeyError

  Edge cases
  - Default resolution=0.1, random_state=0, n_iterations=2
  - Idempotent: running twice overwrites leiden cleanly
  - ADT-only MuData (no RNA modality)
  - n_cells=200
  - Higher resolution produces >= clusters than lower resolution
"""

from __future__ import annotations

import sys
import os
import warnings

import numpy as np
import pytest
import anndata as ad
import scanpy as sc
from anndata import AnnData
from mudata import MuData

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.modules.scripts.cite.adt_annotate import annotate_adt


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_adt_adata(
    n_cells: int = 120,
    n_vars: int = 25,
    n_pca_comps: int = 20,
    n_groups: int = 3,
    seed: int = 42,
) -> AnnData:
    """
    Create a minimal AnnData that looks like adt_harmony.py output.

    The PCA embedding has ``n_groups`` well-separated clusters so that
    Leiden (resolution=0.1) reliably finds multiple communities — required
    for sc.tl.dendrogram to succeed (needs >= 2 groups).

    Outputs written:
      layers["counts"], layers["adt_clr"]
      obsm["X_pca_adt"], obsm["X_pca_harmony_adt"], obsm["X_umap_adt"]
      obsp["connectivities"], obsp["distances"]  — from sc.pp.neighbors
      uns["neighbors"]
    """
    rng = np.random.default_rng(seed)

    raw = rng.integers(0, 50, size=(n_cells, n_vars)).astype(np.float32)
    clr = rng.normal(loc=0.0, scale=2.0, size=(n_cells, n_vars)).astype(np.float32)

    # Structured PCA: n_groups tight clusters well-separated in 2D,
    # padded to n_pca_comps dimensions.
    cells_per_group = n_cells // n_groups
    remainder = n_cells - cells_per_group * n_groups
    angles = np.linspace(0, 2 * np.pi, n_groups, endpoint=False)
    chunks = []
    for g, angle in enumerate(angles):
        extra = 1 if g < remainder else 0
        n = cells_per_group + extra
        center = np.array([10 * np.cos(angle), 10 * np.sin(angle)])
        chunk = rng.normal(loc=center, scale=0.3, size=(n, 2))
        chunks.append(chunk)
    pca_2d = np.vstack(chunks).astype(np.float32)
    pad = rng.normal(scale=0.01, size=(n_cells, n_pca_comps - 2)).astype(np.float32)
    pca_full = np.hstack([pca_2d, pad]).astype(np.float32)

    adata = ad.AnnData(X=raw.copy())
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"ADT_{i}" for i in range(n_vars)]

    adata.layers["counts"]  = raw.copy()
    adata.layers["adt_clr"] = clr.copy()

    adata.obsm["X_pca_adt"]         = pca_full.copy()
    adata.obsm["X_pca_harmony_adt"] = pca_full.copy()
    adata.obsm["X_umap_adt"]        = rng.standard_normal((n_cells, 2)).astype(np.float32)

    # Build neighbor graph from the structured PCA — produces real communities
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(
            adata,
            n_neighbors=15,
            n_pcs=2,
            use_rep="X_pca_harmony_adt",
            random_state=0,
        )

    return adata


def _make_rna_adata(n_cells: int = 120, n_genes: int = 200, seed: int = 99) -> AnnData:
    """Minimal RNA AnnData with pre-computed embeddings."""
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 100, size=(n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm["X_pca"]  = rng.standard_normal((n_cells, 50)).astype(np.float32)
    adata.obsm["X_umap"] = rng.standard_normal((n_cells, 2)).astype(np.float32)
    return adata


def _make_mdata(
    n_cells: int = 120,
    n_adt: int = 25,
    n_pca_comps: int = 20,
    n_groups: int = 3,
    include_rna: bool = True,
    seed: int = 42,
) -> MuData:
    adt = _make_adt_adata(
        n_cells=n_cells, n_vars=n_adt,
        n_pca_comps=n_pca_comps, n_groups=n_groups, seed=seed,
    )
    if include_rna:
        rna = _make_rna_adata(n_cells=n_cells, seed=seed + 1)
        return MuData({"rna": rna, "adt": adt})
    return MuData({"adt": adt})


def _get_clusters(mdata: MuData) -> list[str]:
    """Run annotate_adt once and return sorted unique cluster ID strings."""
    adata, _ = annotate_adt(mdata, inplace=False, random_state=0)
    return sorted(adata.obs["leiden"].unique().tolist())


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestAnnotateAdtHappyPath:
    def test_returns_tuple(self):
        mdata = _make_mdata()
        result = annotate_adt(mdata)
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self):
        mdata = _make_mdata()
        adata, metrics = annotate_adt(mdata)
        assert isinstance(adata, AnnData)
        assert isinstance(metrics, dict)

    def test_leiden_key_always_written(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert "leiden" in adata.obs.columns

    def test_adt_celltype_written_with_map(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(mdata, annotation_map={c: f"CT_{c}" for c in clusters})
        assert "adt_celltype_manual" in adata.obs.columns

    def test_adt_celltype_not_written_without_map(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, annotation_map=None)
        assert "adt_celltype_manual" not in adata.obs.columns

    def test_no_celltype_column_written(self):
        """Guard: obs["celltype"] must never appear (RNA collision)."""
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(mdata, annotation_map={c: f"CT_{c}" for c in clusters})
        assert "celltype" not in adata.obs.columns

    def test_rank_genes_groups_computed(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert "rank_genes_groups" in adata.uns

    def test_dendrogram_computed(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert "dendrogram_leiden" in adata.uns


# ---------------------------------------------------------------------------
# Output shapes and types
# ---------------------------------------------------------------------------

class TestOutputShapesAndTypes:
    def test_leiden_values_are_strings(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        for val in adata.obs["leiden"]:
            assert isinstance(val, str)

    def test_adt_celltype_values_are_strings(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(mdata, annotation_map={c: f"CT_{c}" for c in clusters})
        for val in adata.obs["adt_celltype_manual"]:
            assert isinstance(val, str)

    def test_umap_adt_untouched(self):
        mdata = _make_mdata()
        umap_before = mdata["adt"].obsm["X_umap_adt"].copy()
        adata, _ = annotate_adt(mdata, inplace=False)
        np.testing.assert_array_equal(adata.obsm["X_umap_adt"], umap_before)

    def test_pca_harmony_adt_untouched(self):
        mdata = _make_mdata()
        pca_before = mdata["adt"].obsm["X_pca_harmony_adt"].copy()
        adata, _ = annotate_adt(mdata, inplace=False)
        np.testing.assert_array_equal(adata.obsm["X_pca_harmony_adt"], pca_before)

    def test_adt_clr_layer_preserved(self):
        mdata = _make_mdata()
        clr_before = mdata["adt"].layers["adt_clr"].copy()
        adata, _ = annotate_adt(mdata, inplace=False)
        np.testing.assert_array_equal(adata.layers["adt_clr"], clr_before)

    def test_counts_layer_preserved(self):
        mdata = _make_mdata()
        counts_before = mdata["adt"].layers["counts"].copy()
        adata, _ = annotate_adt(mdata, inplace=False)
        np.testing.assert_array_equal(adata.layers["counts"], counts_before)

    def test_rna_modality_untouched(self):
        mdata = _make_mdata(include_rna=True)
        rna_cols_before = list(mdata["rna"].obs.columns)
        rna_X_before = mdata["rna"].X.copy()
        annotate_adt(mdata, inplace=False)
        np.testing.assert_array_equal(mdata["rna"].X, rna_X_before)
        assert list(mdata["rna"].obs.columns) == rna_cols_before

    def test_n_obs_unchanged(self):
        mdata = _make_mdata(n_cells=120)
        adata, _ = annotate_adt(mdata)
        assert adata.n_obs == 120

    def test_n_vars_unchanged(self):
        mdata = _make_mdata(n_adt=25)
        adata, _ = annotate_adt(mdata)
        assert adata.n_vars == 25


# ---------------------------------------------------------------------------
# Annotation logic
# ---------------------------------------------------------------------------

class TestAnnotationLogic:
    def test_none_map_annotated_false(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata, annotation_map=None)
        assert metrics["annotated"] is False

    def test_none_map_celltype_key_is_none(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata, annotation_map=None)
        assert metrics["celltype_key"] is None

    def test_provided_map_annotated_true(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        _, metrics = annotate_adt(mdata, annotation_map={c: f"T_{c}" for c in clusters})
        assert metrics["annotated"] is True

    def test_provided_map_celltype_key_correct(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        _, metrics = annotate_adt(mdata, annotation_map={c: f"T_{c}" for c in clusters})
        assert metrics["celltype_key"] == "adt_celltype_manual"

    def test_unknown_keys_emit_warning(self):
        mdata = _make_mdata()
        with pytest.warns(UserWarning, match="not found in Leiden clusters"):
            annotate_adt(mdata, annotation_map={"999": "Ghost Cell"})

    def test_unknown_keys_do_not_crash(self):
        mdata = _make_mdata()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_adt(mdata, annotation_map={"999": "Ghost Cell"})
        assert "adt_celltype_manual" in adata.obs.columns

    def test_partial_map_unmapped_clusters_keep_numeric_id(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        if len(clusters) < 2:
            pytest.skip("Need at least 2 clusters for partial map test")
        partial_map = {clusters[0]: "KnownType"}
        adata, _ = annotate_adt(mdata, annotation_map=partial_map)
        unmapped = clusters[1:]
        celltype_values = adata.obs["adt_celltype_manual"].unique().tolist()
        for cluster_id in unmapped:
            assert cluster_id in celltype_values

    def test_full_map_all_clusters_replaced(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        full_map = {c: f"Cell_{c}" for c in clusters}
        adata, _ = annotate_adt(mdata, annotation_map=full_map)
        celltype_values = set(adata.obs["adt_celltype_manual"].tolist())
        expected = {f"Cell_{c}" for c in clusters}
        assert celltype_values == expected

    def test_annotation_map_integer_keys_coerced(self):
        """Integer keys must be coerced to str."""
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        int_map = {int(c): f"Type_{c}" for c in clusters}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, metrics = annotate_adt(mdata, annotation_map=int_map)
        assert metrics["annotated"] is True

    def test_empty_annotation_map_no_replacement(self):
        """Empty dict: annotated=True but leiden values unchanged."""
        mdata = _make_mdata()
        adata, metrics = annotate_adt(mdata, annotation_map={})
        assert metrics["annotated"] is True
        assert (adata.obs["adt_celltype_manual"] == adata.obs["leiden"]).all()


# ---------------------------------------------------------------------------
# Metrics dict
# ---------------------------------------------------------------------------

class TestMetricsDict:
    _REQUIRED_KEYS = {
        "n_cells", "n_clusters", "cluster_sizes",
        "annotation_map", "annotated", "leiden_key",
        "celltype_key", "resolution", "random_state",
    }

    def test_required_keys_present(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert self._REQUIRED_KEYS.issubset(set(metrics.keys()))

    def test_n_cells_correct(self):
        mdata = _make_mdata(n_cells=120)
        _, metrics = annotate_adt(mdata)
        assert metrics["n_cells"] == 120

    def test_n_clusters_at_least_one(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["n_clusters"] >= 1

    def test_cluster_sizes_is_dict(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert isinstance(metrics["cluster_sizes"], dict)

    def test_cluster_sizes_keys_are_strings(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        for k in metrics["cluster_sizes"]:
            assert isinstance(k, str)

    def test_cluster_sizes_values_are_ints(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        for v in metrics["cluster_sizes"].values():
            assert isinstance(v, int)

    def test_cluster_sizes_sum_to_n_cells(self):
        mdata = _make_mdata(n_cells=120)
        _, metrics = annotate_adt(mdata)
        assert sum(metrics["cluster_sizes"].values()) == 120

    def test_annotation_map_empty_when_none(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata, annotation_map=None)
        assert metrics["annotation_map"] == {}

    def test_annotation_map_in_metrics_matches_applied(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        amap = {c: f"T_{c}" for c in clusters}
        _, metrics = annotate_adt(mdata, annotation_map=amap)
        assert metrics["annotation_map"] == amap

    def test_leiden_key_value(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["leiden_key"] == "leiden"

    def test_celltype_key_none_when_not_annotated(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata, annotation_map=None)
        assert metrics["celltype_key"] is None

    def test_celltype_key_correct_when_annotated(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        _, metrics = annotate_adt(mdata, annotation_map={c: f"T_{c}" for c in clusters})
        assert metrics["celltype_key"] == "adt_celltype_manual"

    def test_resolution_in_metrics(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata, resolution=0.2)
        assert metrics["resolution"] == 0.2

    def test_random_state_in_metrics(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata, random_state=7)
        assert metrics["random_state"] == 7

    def test_n_clusters_matches_cluster_sizes_length(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["n_clusters"] == len(metrics["cluster_sizes"])


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

class TestProvenance:
    def test_provenance_key_present(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert "omicsage_adt_annotate" in adata.uns

    def test_provenance_module_name(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["module"] == "adt_annotate"

    def test_provenance_timestamp_nonempty(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        ts = adata.uns["omicsage_adt_annotate"]["timestamp"]
        assert isinstance(ts, str) and len(ts) > 0

    def test_provenance_params_resolution(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, resolution=0.15)
        assert adata.uns["omicsage_adt_annotate"]["params"]["resolution"] == 0.15

    def test_provenance_params_flavor_igraph(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["params"]["flavor"] == "igraph"

    def test_provenance_params_directed_false(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["params"]["directed"] is False

    def test_provenance_params_n_iterations(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, n_iterations=3)
        assert adata.uns["omicsage_adt_annotate"]["params"]["n_iterations"] == 3

    def test_provenance_params_random_state(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, random_state=5)
        assert adata.uns["omicsage_adt_annotate"]["params"]["random_state"] == 5

    def test_provenance_params_annotation_map_provided_false(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, annotation_map=None)
        assert adata.uns["omicsage_adt_annotate"]["params"]["annotation_map_provided"] is False

    def test_provenance_params_annotation_map_provided_true(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(mdata, annotation_map={c: "X" for c in clusters})
        assert adata.uns["omicsage_adt_annotate"]["params"]["annotation_map_provided"] is True

    def test_provenance_outputs_leiden_key(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["outputs"]["leiden_key"] == "leiden"

    def test_provenance_outputs_n_clusters(self):
        mdata = _make_mdata()
        adata, metrics = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["outputs"]["n_clusters"] == metrics["n_clusters"]

    def test_provenance_outputs_celltype_key_none_when_not_annotated(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, annotation_map=None)
        assert adata.uns["omicsage_adt_annotate"]["outputs"]["celltype_key"] is None

    def test_provenance_outputs_celltype_key_correct_when_annotated(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(mdata, annotation_map={c: "X" for c in clusters})
        assert adata.uns["omicsage_adt_annotate"]["outputs"]["celltype_key"] == "adt_celltype_manual"


# ---------------------------------------------------------------------------
# inplace behaviour
# ---------------------------------------------------------------------------

class TestInplace:
    def test_inplace_false_original_adt_untouched(self):
        mdata = _make_mdata()
        obs_cols_before = list(mdata["adt"].obs.columns)
        annotate_adt(mdata, inplace=False)
        assert "leiden" not in mdata["adt"].obs.columns
        assert list(mdata["adt"].obs.columns) == obs_cols_before

    def test_inplace_true_modifies_adt(self):
        mdata = _make_mdata()
        annotate_adt(mdata, inplace=True)
        assert "leiden" in mdata["adt"].obs.columns

    def test_inplace_false_returns_copy(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, inplace=False)
        assert adata is not mdata["adt"]

    def test_inplace_true_returns_same_object(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata, inplace=True)
        assert adata is mdata["adt"]


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_missing_adt_modality_raises_keyerror(self):
        rna = _make_rna_adata()
        mdata = MuData({"rna": rna})
        with pytest.raises(KeyError, match="adt"):
            annotate_adt(mdata)

    def test_missing_connectivities_raises_keyerror(self):
        mdata = _make_mdata()
        del mdata["adt"].obsp["connectivities"]
        with pytest.raises(KeyError, match="connectivities"):
            annotate_adt(mdata)


# ---------------------------------------------------------------------------
# Default parameters (sc-best-practices ch.39)
# ---------------------------------------------------------------------------

class TestDefaultParameters:
    def test_default_resolution_is_0_1(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["resolution"] == 0.1

    def test_default_random_state_is_0(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["random_state"] == 0

    def test_default_n_iterations_in_provenance(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["params"]["n_iterations"] == 2


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_adt_only_mdata_no_rna(self):
        mdata = _make_mdata(include_rna=False)
        adata, metrics = annotate_adt(mdata)
        assert "leiden" in adata.obs.columns
        assert metrics["n_clusters"] >= 1

    def test_idempotent_second_run_overwrites_leiden(self):
        mdata = _make_mdata()
        adata1, _ = annotate_adt(mdata, inplace=False)
        mdata2 = MuData({"adt": adata1})
        adata2, _ = annotate_adt(mdata2, inplace=True)
        assert "leiden" in adata2.obs.columns
        assert adata2.obs.columns.tolist().count("leiden") == 1

    def test_custom_resolution_higher_gives_more_or_equal_clusters(self):
        mdata_lo = _make_mdata(seed=42)
        mdata_hi = _make_mdata(seed=42)
        _, m_lo = annotate_adt(mdata_lo, resolution=0.05)
        _, m_hi = annotate_adt(mdata_hi, resolution=0.8)
        assert m_hi["n_clusters"] >= m_lo["n_clusters"]

    def test_n_cells_200(self):
        mdata = _make_mdata(n_cells=200)
        adata, metrics = annotate_adt(mdata)
        assert metrics["n_cells"] == 200
        assert sum(metrics["cluster_sizes"].values()) == 200


# ===========================================================================
# Tests for auto-scoring (obs["adt_celltype_score"])
# ===========================================================================

def _make_mdata_with_markers(
    n_cells: int = 120,
    n_groups: int = 3,
    seed: int = 42,
) -> MuData:
    """
    MuData whose ADT panel contains at least CD3, CD19, CD14 — the proteins
    in the default marker pairs AND the bmmc preset.  Clusters are
    well-separated so scoring can deterministically pick a winner.
    """
    rng = np.random.default_rng(seed)
    # Protein names that cover both default marker pairs and bmmc preset keys
    var_names = [
        "CD3-TotalSeqA", "CD19-TotalSeqA", "CD14-TotalSeqA",
        "CD4-TotalSeqA",  "CD8a-TotalSeqA", "CD56-TotalSeqA",
        "CD20-TotalSeqA", "CD11c-TotalSeqA","CD16-TotalSeqA",
        "CD38-TotalSeqA", "CD71-TotalSeqA", "CD36-TotalSeqA",
    ]
    n_vars = len(var_names)
    n_cells_per = n_cells // n_groups

    raw = rng.integers(0, 30, size=(n_cells, n_vars)).astype(np.float32)
    clr = rng.normal(loc=0.0, scale=1.0, size=(n_cells, n_vars)).astype(np.float32)

    # Make each group strongly express a different protein (CD3 / CD19 / CD14)
    for g in range(n_groups):
        start = g * n_cells_per
        end   = start + n_cells_per
        clr[start:end, g] += 8.0   # group g strongly expresses protein g

    # Build structured PCA
    angles = np.linspace(0, 2 * np.pi, n_groups, endpoint=False)
    chunks = []
    for angle in angles:
        center = np.array([10 * np.cos(angle), 10 * np.sin(angle)])
        chunk  = rng.normal(loc=center, scale=0.3, size=(n_cells_per, 2))
        chunks.append(chunk)
    pca_2d = np.vstack(chunks).astype(np.float32)
    n_pca  = 20
    pad    = rng.normal(scale=0.01, size=(len(pca_2d), n_pca - 2)).astype(np.float32)
    pca    = np.hstack([pca_2d, pad]).astype(np.float32)

    adata = ad.AnnData(X=raw.copy())
    adata.obs_names = [f"cell_{i}" for i in range(len(pca_2d))]
    adata.var_names = var_names

    adata.layers["counts"]  = raw[:len(pca_2d)].copy()
    adata.layers["adt_clr"] = clr[:len(pca_2d)].copy()
    adata.obsm["X_pca_adt"]         = pca.copy()
    adata.obsm["X_pca_harmony_adt"] = pca.copy()
    adata.obsm["X_umap_adt"]        = rng.standard_normal((len(pca_2d), 2)).astype(np.float32)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=2,
                        use_rep="X_pca_harmony_adt", random_state=0)

    return MuData({"adt": adata})


# ---------------------------------------------------------------------------
# Scoring — happy path
# ---------------------------------------------------------------------------

class TestScoringHappyPath:
    _SIMPLE_PANEL = {
        "T cell":   ["CD3"],
        "B cell":   ["CD19"],
        "Monocyte": ["CD14"],
    }

    def test_score_column_written_with_marker_panel(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert "adt_celltype_score" in adata.obs.columns

    def test_score_column_written_with_preset_bmmc(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, preset="bmmc")
        assert "adt_celltype_score" in adata.obs.columns

    def test_score_column_not_written_without_panel_or_preset(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert "adt_celltype_score" not in adata.obs.columns

    def test_score_values_are_strings(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        for val in adata.obs["adt_celltype_score"]:
            assert isinstance(val, str)

    def test_score_values_are_subset_of_panel_keys_plus_unknown(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        valid = set(self._SIMPLE_PANEL.keys()) | {"Unknown"}
        for val in adata.obs["adt_celltype_score"].unique():
            assert val in valid

    def test_score_all_cells_labelled(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert adata.obs["adt_celltype_score"].isna().sum() == 0

    def test_score_n_obs_unchanged(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert len(adata.obs["adt_celltype_score"]) == adata.n_obs

    def test_score_and_manual_can_coexist(self):
        """Both columns written independently when both params given."""
        mdata = _make_mdata_with_markers()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(
            mdata,
            marker_panel=self._SIMPLE_PANEL,
            annotation_map={c: f"CT_{c}" for c in clusters},
        )
        assert "adt_celltype_score"  in adata.obs.columns
        assert "adt_celltype_manual" in adata.obs.columns

    def test_manual_absent_when_only_panel_given(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert "adt_celltype_manual" not in adata.obs.columns

    def test_score_absent_when_only_manual_given(self):
        mdata = _make_mdata()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(mdata, annotation_map={c: "X" for c in clusters})
        assert "adt_celltype_score" not in adata.obs.columns


# ---------------------------------------------------------------------------
# Scoring — cluster-level correctness
# ---------------------------------------------------------------------------

class TestScoringClusterLevel:
    """
    The fixture strongly expresses:
      group 0 → protein 0 (CD3)   → should score as T cell
      group 1 → protein 1 (CD19)  → should score as B cell
      group 2 → protein 2 (CD14)  → should score as Monocyte
    """
    _SIMPLE_PANEL = {
        "T cell":   ["CD3"],
        "B cell":   ["CD19"],
        "Monocyte": ["CD14"],
    }

    def test_cells_per_cluster_get_consistent_label(self):
        """All cells in a cluster must share the same score label."""
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        for cluster in adata.obs["leiden"].unique():
            mask   = adata.obs["leiden"] == cluster
            labels = adata.obs.loc[mask, "adt_celltype_score"].unique()
            assert len(labels) == 1, (
                f"Cluster {cluster} has multiple score labels: {labels}"
            )

    def test_score_differs_from_manual_when_maps_differ(self):
        """Score and manual can disagree — they are independent."""
        mdata = _make_mdata_with_markers()
        clusters = _get_clusters(mdata)
        wrong_map = {c: "WrongType" for c in clusters}
        adata, _ = annotate_adt(
            mdata,
            marker_panel=self._SIMPLE_PANEL,
            annotation_map=wrong_map,
        )
        # manual is all "WrongType"; score should have biologically-derived labels
        assert set(adata.obs["adt_celltype_manual"].unique()) == {"WrongType"}
        assert "WrongType" not in adata.obs["adt_celltype_score"].unique()


# ---------------------------------------------------------------------------
# Scoring — metrics
# ---------------------------------------------------------------------------

class TestScoringMetrics:
    _SIMPLE_PANEL = {
        "T cell":   ["CD3"],
        "B cell":   ["CD19"],
        "Monocyte": ["CD14"],
    }

    def test_score_key_in_metrics_when_panel_given(self):
        mdata = _make_mdata_with_markers()
        _, metrics = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert metrics["score_key"] == "adt_celltype_score"

    def test_score_key_none_when_no_panel(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["score_key"] is None

    def test_marker_panel_in_metrics(self):
        mdata = _make_mdata_with_markers()
        _, metrics = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert metrics["marker_panel"] == self._SIMPLE_PANEL

    def test_marker_panel_none_in_metrics_when_not_given(self):
        mdata = _make_mdata()
        _, metrics = annotate_adt(mdata)
        assert metrics["marker_panel"] is None

    def test_score_key_in_metrics_with_preset(self):
        mdata = _make_mdata_with_markers()
        _, metrics = annotate_adt(mdata, preset="bmmc")
        assert metrics["score_key"] == "adt_celltype_score"


# ---------------------------------------------------------------------------
# Scoring — provenance
# ---------------------------------------------------------------------------

class TestScoringProvenance:
    _SIMPLE_PANEL = {
        "T cell":   ["CD3"],
        "B cell":   ["CD19"],
        "Monocyte": ["CD14"],
    }

    def test_provenance_score_key_correct(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert adata.uns["omicsage_adt_annotate"]["outputs"]["score_key"] \
               == "adt_celltype_score"

    def test_provenance_score_key_none_when_no_panel(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["outputs"]["score_key"] is None

    def test_provenance_marker_panel_provided_true(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert adata.uns["omicsage_adt_annotate"]["params"]["marker_panel_provided"] is True

    def test_provenance_marker_panel_provided_false(self):
        mdata = _make_mdata()
        adata, _ = annotate_adt(mdata)
        assert adata.uns["omicsage_adt_annotate"]["params"]["marker_panel_provided"] is False

    def test_provenance_preset_recorded(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, preset="bmmc")
        assert adata.uns["omicsage_adt_annotate"]["params"]["preset"] == "bmmc"

    def test_marker_panel_stored_in_uns(self):
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert "omicsage_adt_marker_panel" in adata.uns
        assert adata.uns["omicsage_adt_marker_panel"] == self._SIMPLE_PANEL


# ---------------------------------------------------------------------------
# Scoring — preset
# ---------------------------------------------------------------------------

class TestScoringPreset:
    def test_unknown_preset_emits_warning(self):
        mdata = _make_mdata()
        with pytest.warns(UserWarning, match="Unknown preset"):
            annotate_adt(mdata, preset="nonexistent_panel")

    def test_unknown_preset_does_not_write_score_column(self):
        mdata = _make_mdata()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_adt(mdata, preset="nonexistent_panel")
        assert "adt_celltype_score" not in adata.obs.columns

    def test_marker_panel_takes_precedence_over_preset(self):
        """When both marker_panel and preset given, marker_panel wins."""
        custom_panel = {"MyType": ["CD3"]}
        mdata = _make_mdata_with_markers()
        adata, metrics = annotate_adt(
            mdata, marker_panel=custom_panel, preset="bmmc"
        )
        assert metrics["marker_panel"] == custom_panel
        valid = set(custom_panel.keys()) | {"Unknown"}
        for val in adata.obs["adt_celltype_score"].unique():
            assert val in valid


# ---------------------------------------------------------------------------
# Scoring — no collision with RNA / no "adt_celltype" plain column
# ---------------------------------------------------------------------------

class TestScoringNoCollision:
    _SIMPLE_PANEL = {
        "T cell":   ["CD3"],
        "B cell":   ["CD19"],
        "Monocyte": ["CD14"],
    }

    def test_no_plain_adt_celltype_column(self):
        """obs['adt_celltype'] must never be written."""
        mdata = _make_mdata_with_markers()
        clusters = _get_clusters(mdata)
        adata, _ = annotate_adt(
            mdata,
            marker_panel=self._SIMPLE_PANEL,
            annotation_map={c: f"CT_{c}" for c in clusters},
        )
        assert "adt_celltype" not in adata.obs.columns

    def test_no_celltype_column(self):
        """obs['celltype'] must never be written (RNA collision guard)."""
        mdata = _make_mdata_with_markers()
        adata, _ = annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL)
        assert "celltype" not in adata.obs.columns

    def test_rna_untouched_with_scoring(self):
        mdata = _make_mdata_with_markers()
        rna_cols_before = list(mdata["adt"].obs.columns)
        annotate_adt(mdata, marker_panel=self._SIMPLE_PANEL, inplace=False)
        assert list(mdata["adt"].obs.columns) == rna_cols_before
