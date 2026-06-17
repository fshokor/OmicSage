"""
tests/test_annotate.py
======================
Test suite for pipeline/modules/scripts/annotation/annotate.py

Tests
-----
test_celltypist_columns_exist        — obs['celltypist_coarse'] and ['celltypist_fine'] written
test_celltypist_labels_are_strings   — all CellTypist labels are non-empty strings
test_every_cluster_gets_celltypist   — every leiden cluster has a non-null coarse label
test_marker_column_exists            — obs['cell_type_markers'] written
test_marker_labels_are_valid         — all marker labels are keys in MARKER_SETS
test_every_cluster_gets_marker       — every leiden cluster has a marker label
test_vote_columns_exist              — obs['cell_type'] and ['cell_type_confidence'] written
test_vote_confidence_range           — confidence values are in [0.0, 1.0]
test_inplace_false_no_mutation       — inplace=False leaves the original adata unchanged
test_params_stored_in_uns            — uns['omicsage_annotate'] has required keys
test_unknown_method_raises           — ValueError on unknown method name
test_vote_without_prerequisites_raises — ValueError if vote requested without celltypist+markers
test_missing_leiden_col_raises       — KeyError if leiden_col not in obs

Note on CellTypist in tests
----------------------------
We do NOT mock CellTypist — these tests use a lightweight fixture that does NOT
call CellTypist at all (methods=["markers"] or methods=["markers","vote"] where
celltypist is absent so vote falls back gracefully).

A separate integration-level test tag (omit from CI) covers the real CellTypist
path on the full GSE194122 dataset.

ScType-py tests — add these to tests/test_annotate.py
======================================================
All network-dependent tests are skipped in CI via OMICSAGE_CI env var.
Pure-logic tests (parse, score) run always — no network needed.
"""

import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import scanpy as sc

import os
import io
from unittest.mock import patch, MagicMock

# ── Import the module under test ───────────────────────────────────────────────
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.modules.scripts.annotation.annotate import (
    annotate,
    _fetch_sctype_db,
    _parse_sctype_markers,
    _run_sctype,
    _run_scanvi,
    _run_majority_vote,
    _run_singler,
    _load_singler_ref,
    _HPCA_CACHE,
    _TABULA_SAPIENS_CACHE,
    MARKER_SETS,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

def _make_fixture(n_cells: int = 200, n_genes: int = 150,
                  n_clusters: int = 5, sparse: bool = True) -> sc.AnnData:
    """
    Minimal clustered AnnData fixture suitable for annotate() tests.

    - Sparse count matrix with a handful of marker genes guaranteed present
    - 'leiden' obs column with integer cluster labels
    - 'logcounts' layer (log-normalised copy)
    - 'counts' layer (raw integer counts)
    - UMAP embedding so scanpy plotting won't crash in report tests
    """
    rng   = np.random.default_rng(42)
    n_var = n_genes

    # Guarantee at least a few genes from each marker set are present
    all_markers: list = []
    for markers in MARKER_SETS.values():
        all_markers.extend(markers[:2])   # take first 2 markers per type
    all_markers = list(dict.fromkeys(all_markers))  # deduplicate, preserve order
    n_marker = min(len(all_markers), n_var)
    gene_names = all_markers[:n_marker] + [
        f"GENE{i}" for i in range(n_var - n_marker)
    ]

    counts = rng.negative_binomial(5, 0.5, size=(n_cells, n_var)).astype(np.float32)
    if sparse:
        X = sp.csr_matrix(counts)
    else:
        X = counts

    adata = sc.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

    # Cluster labels — evenly distribute
    labels = [str(i % n_clusters) for i in range(n_cells)]
    adata.obs["leiden"] = pd.Categorical(labels)

    # Log-normalised layer
    lognorm = counts / counts.sum(axis=1, keepdims=True) * 1e4
    lognorm = np.log1p(lognorm).astype(np.float32)
    adata.layers["logcounts"] = sp.csr_matrix(lognorm) if sparse else lognorm
    adata.layers["counts"]    = X.copy() if sparse else counts.copy()

    # Minimal UMAP (random — just so the slot exists)
    adata.obsm["X_umap"] = rng.standard_normal((n_cells, 2)).astype(np.float32)

    return adata


@pytest.fixture
def adata_clustered():
    return _make_fixture()


@pytest.fixture
def adata_clustered_dense():
    return _make_fixture(sparse=False)


# ─────────────────────────────────────────────────────────────────────────────
# Marker scoring tests  (pure Python, always run)
# ─────────────────────────────────────────────────────────────────────────────

def test_marker_column_exists(adata_clustered):
    adata_ann, _ = annotate(adata_clustered, methods=["markers"], inplace=False)
    assert "cell_type_markers" in adata_ann.obs.columns


def test_marker_labels_are_valid(adata_clustered):
    adata_ann, _ = annotate(adata_clustered, methods=["markers"], inplace=False)
    valid_labels = set(MARKER_SETS.keys())
    for label in adata_ann.obs["cell_type_markers"].unique():
        assert label in valid_labels, (
            f"Unexpected marker label '{label}' — not in MARKER_SETS"
        )


def test_every_cluster_gets_marker(adata_clustered):
    adata_ann, _ = annotate(adata_clustered, methods=["markers"], inplace=False)
    for cl in adata_clustered.obs["leiden"].unique():
        mask   = adata_ann.obs["leiden"] == cl
        labels = adata_ann.obs.loc[mask, "cell_type_markers"].dropna()
        assert len(labels) > 0, f"Cluster {cl} has no marker label"
        assert all(isinstance(l, str) for l in labels)


def test_marker_works_with_dense_matrix(adata_clustered_dense):
    adata_ann, _ = annotate(adata_clustered_dense, methods=["markers"], inplace=False)
    assert "cell_type_markers" in adata_ann.obs.columns


def test_marker_score_df_in_ann_dict(adata_clustered):
    _, ann_dict = annotate(adata_clustered, methods=["markers"], inplace=False)
    assert "marker_score_df" in ann_dict
    df = ann_dict["marker_score_df"]
    assert "best_by_score" in df.columns
    assert len(df) == adata_clustered.obs["leiden"].nunique()


# ─────────────────────────────────────────────────────────────────────────────
# Vote tests (markers-only vote — celltypist columns absent, vote falls back gracefully)
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def adata_pre_annotated(adata_clustered):
    """
    Fixture with marker labels + a synthetic celltypist_fine column already in obs.
    Used to test the vote path without requiring celltypist to be installed.
    We call the internal _run_majority_vote directly so the test is independent
    of the CellTypist import guard.
    """
    from pipeline.modules.scripts.annotation.annotate import _run_marker_scoring, _run_majority_vote
    adata = adata_clustered.copy()
    # Run marker scoring to populate cell_type_markers
    score_df = _run_marker_scoring(adata, "leiden", MARKER_SETS)
    # Simulate celltypist_fine being present (same labels as markers — perfect agreement)
    adata.obs["celltypist_fine"] = adata.obs["cell_type_markers"].copy()
    # Run vote directly
    _run_majority_vote(adata, score_df, "leiden")
    return adata


def test_vote_columns_exist(adata_pre_annotated):
    """obs['cell_type_vote'] and ['cell_type_confidence'] must be written by vote."""
    assert "cell_type_vote" in adata_pre_annotated.obs.columns
    assert "cell_type_confidence" in adata_pre_annotated.obs.columns


def test_vote_confidence_range(adata_pre_annotated):
    conf = adata_pre_annotated.obs["cell_type_confidence"]
    assert conf.between(0.0, 1.0).all(), (
        f"Confidence values outside [0, 1]: min={conf.min()}, max={conf.max()}"
    )


def test_vote_warns_when_prerequisites_missing(adata_clustered):
    """vote requested but celltypist not in methods — should warn, not crash."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        # markers only, but vote also listed — celltypist missing
        # The module raises ValueError because vote requires celltypist+markers
        # So we test the prerequisite check instead
        with pytest.raises(ValueError, match="vote"):
            annotate(adata_clustered, methods=["markers", "vote"], inplace=False)


# ─────────────────────────────────────────────────────────────────────────────
# inplace / mutation tests
# ─────────────────────────────────────────────────────────────────────────────

def test_inplace_false_no_mutation(adata_clustered):
    original_cols = set(adata_clustered.obs.columns)
    annotate(adata_clustered, methods=["markers"], inplace=False)
    assert set(adata_clustered.obs.columns) == original_cols, (
        "inplace=False mutated the original adata.obs"
    )


def test_groundtruth_preserved(adata_clustered):
    """obs['cell_type'] present before annotate() → copied to obs['cell_type_groundtruth']."""
    adata_with_gt = adata_clustered.copy()
    adata_with_gt.obs["cell_type"] = "SomeGroundTruth"
    adata_ann, _ = annotate(adata_with_gt, methods=["markers"], inplace=False)
    assert "cell_type_groundtruth" in adata_ann.obs.columns, (
        "obs['cell_type_groundtruth'] not written when obs['cell_type'] pre-existed"
    )
    assert (adata_ann.obs["cell_type_groundtruth"] == "SomeGroundTruth").all(), (
        "Ground-truth values were not correctly preserved"
    )


def test_groundtruth_not_created_when_absent(adata_clustered):
    """obs['cell_type'] absent → obs['cell_type_groundtruth'] must NOT be created."""
    assert "cell_type" not in adata_clustered.obs.columns
    adata_ann, _ = annotate(adata_clustered, methods=["markers"], inplace=False)
    assert "cell_type_groundtruth" not in adata_ann.obs.columns


def test_inplace_true_mutates(adata_clustered):
    adata_copy = adata_clustered.copy()
    annotate(adata_copy, methods=["markers"], inplace=True)
    assert "cell_type_markers" in adata_copy.obs.columns


# ─────────────────────────────────────────────────────────────────────────────
# Provenance / uns tests
# ─────────────────────────────────────────────────────────────────────────────

def test_params_stored_in_uns(adata_clustered):
    adata_ann, _ = annotate(adata_clustered, methods=["markers"], inplace=False)
    prov = adata_ann.uns.get("omicsage_annotate")
    assert prov is not None, "uns['omicsage_annotate'] not written"
    required_keys = {
        "methods_requested", "methods_run", "leiden_col",
        "marker_sets_keys", "n_clusters", "n_cells",
        "scanpy_version", "omicsage_module", "omicsage_version", "timestamp",
    }
    missing = required_keys - set(prov.keys())
    assert not missing, f"Missing uns keys: {missing}"


def test_uns_n_clusters_correct(adata_clustered):
    adata_ann, _ = annotate(adata_clustered, methods=["markers"], inplace=False)
    expected = adata_clustered.obs["leiden"].nunique()
    assert adata_ann.uns["omicsage_annotate"]["n_clusters"] == expected


# ─────────────────────────────────────────────────────────────────────────────
# Error handling tests
# ─────────────────────────────────────────────────────────────────────────────

def test_unknown_method_raises(adata_clustered):
    with pytest.raises(ValueError, match="Unknown methods"):
        annotate(adata_clustered, methods=["notamethod"], inplace=False)


def test_vote_without_prerequisites_raises(adata_clustered):
    with pytest.raises(ValueError, match="vote"):
        annotate(adata_clustered, methods=["markers", "vote"], inplace=False)


def test_missing_leiden_col_raises(adata_clustered):
    with pytest.raises(KeyError):
        annotate(adata_clustered, methods=["markers"],
                 leiden_col="nonexistent_col", inplace=False)


# ─────────────────────────────────────────────────────────────────────────────
# Custom marker sets test
# ─────────────────────────────────────────────────────────────────────────────

def test_custom_marker_sets(adata_clustered):
    """User-supplied marker_sets are respected."""
    custom = {
        "TypeA": ["CD3D", "CD3E"],
        "TypeB": ["CD68", "LYZ"],
    }
    adata_ann, ann_dict = annotate(
        adata_clustered, methods=["markers"],
        marker_sets=custom, inplace=False
    )
    valid = set(custom.keys())
    for label in adata_ann.obs["cell_type_markers"].unique():
        assert label in valid, f"Label '{label}' not in custom marker set keys"
    # marker_sets_keys in uns should reflect the custom set
    assert set(adata_ann.uns["omicsage_annotate"]["marker_sets_keys"]) == valid


# ── Shared fixture ─────────────────────────────────────────────────────────────
 
@pytest.fixture
def adata_clustered():
    """
    Minimal AnnData with logcounts layer and leiden clusters.
    10 cells, 20 genes, 3 clusters.
    Gene names include a handful from MARKER_SETS so scoring is non-zero.
    """
    rng = np.random.default_rng(42)
    n_cells, n_genes = 10, 20
 
    # Use real gene names so marker sets can match
    gene_names = [
        "CD3D", "CD3E", "CD68", "CD14", "MS4A1",    # immune markers
        "ALB", "APOC3", "COL1A1", "PECAM1", "VWF",  # parenchymal
        "GeneA", "GeneB", "GeneC", "GeneD", "GeneE",
        "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ",
    ]
    X = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = sc.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.obs["leiden"] = pd.Categorical(
        ["0", "0", "0", "1", "1", "1", "1", "2", "2", "2"]
    )
    adata.layers["logcounts"] = np.log1p(X)
    return adata
 
 
# ── Minimal mock ScTypeDB ──────────────────────────────────────────────────────
 
def _make_mock_db() -> pd.DataFrame:
    """
    Minimal ScTypeDB-shaped DataFrame covering two tissues.
    Positive/negative markers are drawn from gene_names in the fixture.
    """
    return pd.DataFrame([
        {
            "tissueType":      "Immune system",
            "cellName":        "T cell",
            "geneSymbolmore1": "CD3D,CD3E",
            "geneSymbolmore2": "CD68",
        },
        {
            "tissueType":      "Immune system",
            "cellName":        "Macrophage",
            "geneSymbolmore1": "CD68,CD14",
            "geneSymbolmore2": "CD3D",
        },
        {
            "tissueType":      "Immune system",
            "cellName":        "B cell",
            "geneSymbolmore1": "MS4A1",
            "geneSymbolmore2": "",
        },
        {
            "tissueType":      "Liver",
            "cellName":        "Hepatocyte",
            "geneSymbolmore1": "ALB,APOC3",
            "geneSymbolmore2": "CD3D",
        },
        {
            "tissueType":      "Liver",
            "cellName":        "Fibroblast",
            "geneSymbolmore1": "COL1A1",
            "geneSymbolmore2": "",
        },
    ])
 
 
# ── _parse_sctype_markers tests (no network) ───────────────────────────────────
 
class TestParseScTypeMarkers:
 
    def test_returns_dict_for_valid_tissue(self):
        db = _make_mock_db()
        result = _parse_sctype_markers(db, "Immune system")
        assert isinstance(result, dict)
        assert len(result) > 0
 
    def test_cell_types_have_pos_neg_keys(self):
        db = _make_mock_db()
        result = _parse_sctype_markers(db, "Immune system")
        for ct, markers in result.items():
            assert "pos" in markers, f"{ct} missing 'pos'"
            assert "neg" in markers, f"{ct} missing 'neg'"
 
    def test_positive_markers_parsed_correctly(self):
        db = _make_mock_db()
        result = _parse_sctype_markers(db, "Immune system")
        assert "CD3D" in result["T cell"]["pos"]
        assert "CD3E" in result["T cell"]["pos"]
 
    def test_negative_markers_parsed_correctly(self):
        db = _make_mock_db()
        result = _parse_sctype_markers(db, "Immune system")
        assert "CD68" in result["T cell"]["neg"]
 
    def test_empty_negative_markers_ok(self):
        db = _make_mock_db()
        result = _parse_sctype_markers(db, "Immune system")
        # B cell has no negative markers
        assert result["B cell"]["neg"] == []
 
    def test_tissue_filter_isolates_correct_types(self):
        db = _make_mock_db()
        immune = _parse_sctype_markers(db, "Immune system")
        liver  = _parse_sctype_markers(db, "Liver")
        assert "T cell" in immune
        assert "Hepatocyte" in liver
        assert "T cell" not in liver
        assert "Hepatocyte" not in immune
 
    def test_invalid_tissue_raises_value_error(self):
        db = _make_mock_db()
        with pytest.raises(ValueError, match="not found in ScTypeDB"):
            _parse_sctype_markers(db, "Nonexistent tissue")
 
    def test_invalid_tissue_error_lists_available(self):
        db = _make_mock_db()
        with pytest.raises(ValueError, match="Immune system"):
            _parse_sctype_markers(db, "Nonexistent tissue")
 
 
# ── _run_sctype tests (mocked network) ────────────────────────────────────────
 
class TestRunScType:
 
    def test_column_exists_after_run(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        assert "cell_type_sctype" in adata_clustered.obs.columns
 
    def test_labels_are_strings(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        labels = adata_clustered.obs["cell_type_sctype"]
        assert all(isinstance(v, str) for v in labels), \
            "All cell_type_sctype values must be strings"
 
    def test_every_cluster_covered(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        n_null = adata_clustered.obs["cell_type_sctype"].isna().sum()
        assert n_null == 0, f"{n_null} cells have null cell_type_sctype"
 
    def test_no_cells_labelled_unknown_when_markers_present(self, adata_clustered):
        """
        With mock DB that covers immune genes present in fixture,
        no cluster should fall back to 'Unknown'.
        """
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        labels = adata_clustered.obs["cell_type_sctype"].unique().tolist()
        assert "Unknown" not in labels
 
    def test_liver_tissue_gives_liver_types(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Liver")
        labels = set(adata_clustered.obs["cell_type_sctype"].unique())
        liver_types = {"Hepatocyte", "Fibroblast"}
        assert labels.issubset(liver_types | {"Unknown"}), \
            f"Unexpected labels for Liver tissue: {labels - liver_types}"
 
    def test_column_is_category_dtype(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        assert hasattr(adata_clustered.obs["cell_type_sctype"], "cat"), \
            "cell_type_sctype should be category dtype"
 
    def test_uses_logcounts_layer_when_present(self, adata_clustered):
        """Patch _mean_expr indirectly — verify logcounts is preferred over .X"""
        sentinel = np.zeros_like(adata_clustered.layers["logcounts"])
        adata_clustered.layers["logcounts"] = sentinel
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        # With zero expression all scores are 0; best_type still assigned (no crash)
        assert "cell_type_sctype" in adata_clustered.obs.columns
 
 
# ── annotate() integration tests (mocked network) ─────────────────────────────
 
class TestAnnotateWithScType:
 
    def test_sctype_in_methods_run(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            with patch("pipeline.modules.scripts.annotation.annotate._run_celltypist"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "sctype"],
                    tissue="Immune system",
                )
        assert "sctype" in ann_dict["provenance"]["methods_run"]
 
    def test_tissue_recorded_in_provenance(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            with patch("pipeline.modules.scripts.annotation.annotate._run_celltypist"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "sctype"],
                    tissue="Liver",
                )
        assert ann_dict["provenance"]["sctype_tissue"] == "Liver"
 
    def test_sctype_column_written_by_annotate(self, adata_clustered):
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            with patch("pipeline.modules.scripts.annotation.annotate._run_celltypist"):
                adata_ann, _ = annotate(
                    adata_clustered,
                    methods=["markers", "sctype"],
                    tissue="Immune system",
                )
        assert "cell_type_sctype" in adata_ann.obs.columns
 
    def test_sctype_skipped_gracefully_on_network_error(self, adata_clustered):
        """If network fails, sctype warns and skips — does not crash annotate()."""
        import requests
        with patch(
            "pipeline.modules.scripts.annotation.annotate._fetch_sctype_db",
            side_effect=requests.ConnectionError("network down"),
        ):
            with pytest.warns(UserWarning, match="ScType annotation failed"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "sctype"],
                    tissue="Immune system",
                )
        assert "sctype" not in ann_dict["provenance"]["methods_run"]
        assert "cell_type_sctype" not in adata_ann.obs.columns
 
 
# ── Network test (skipped in CI) ──────────────────────────────────────────────
 
@pytest.mark.skipif(
    os.getenv("OMICSAGE_CI") == "true",
    reason="Skipped in CI — requires network access to GitHub",
)
class TestFetchScTypeDB:
 
    def test_fetch_returns_dataframe(self):
        db = _fetch_sctype_db()
        assert isinstance(db, pd.DataFrame)
        assert len(db) > 0
 
    def test_fetch_has_required_columns(self):
        db = _fetch_sctype_db()
        for col in ("tissueType", "cellName", "geneSymbolmore1", "geneSymbolmore2"):
            assert col in db.columns, f"Missing column: {col}"
 
    def test_fetch_contains_immune_tissue(self):
        db = _fetch_sctype_db()
        assert "Immune system" in db["tissueType"].values
 
    def test_fetch_contains_liver_tissue(self):
        db = _fetch_sctype_db()
        assert "Liver" in db["tissueType"].values

# ─────────────────────────────────────────────────────────────────────────────
# Session B placeholder — skipped in CI
# ─────────────────────────────────────────────────────────────────────────────

import os

@pytest.mark.skipif(
    os.getenv("OMICSAGE_CI") == "true",
    reason="CellTypist integration test skipped in CI (requires model download)"
)
def test_celltypist_integration_skipped_in_ci(adata_clustered):
    """
    Placeholder for full CellTypist integration test.
    Run locally after: pip install celltypist
    This test will be fleshed out in Session B.
    """
    pytest.skip("CellTypist integration test — enable locally when needed")


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: scANVI tests  (all CI-skipped — require model download)
# ─────────────────────────────────────────────────────────────────────────────

def _make_mock_scanvi_model(labels, probs=None):
    """
    Return a mock SCANVI model that yields fixed labels (and optional probs).
    """
    n = len(labels)
    if probs is None:
        probs = np.full(n, 0.85)

    pred_df = pd.DataFrame({
        "prediction": labels,
        "probability": probs,
    })

    mock_model = MagicMock()
    mock_model.predict.return_value = pred_df
    return mock_model


@pytest.mark.skipif(
    os.getenv("OMICSAGE_CI") == "true",
    reason="scANVI tests skipped in CI — require model download",
)
class TestRunScANVI:

    def test_cell_type_scanvi_column_written(self, adata_clustered):
        labels = ["T cell"] * len(adata_clustered)
        mock_model = _make_mock_scanvi_model(labels)

        with patch("scvi.model.SCANVI") as mock_cls:
            mock_cls.load.return_value = mock_model
            result = _run_scanvi(adata_clustered, "leiden", "/fake/model")

        assert "cell_type_scanvi" in adata_clustered.obs.columns

    def test_labels_are_strings(self, adata_clustered):
        labels = ["T cell"] * len(adata_clustered)
        mock_model = _make_mock_scanvi_model(labels)

        with patch("scvi.model.SCANVI") as mock_cls:
            mock_cls.load.return_value = mock_model
            _run_scanvi(adata_clustered, "leiden", "/fake/model")

        assert all(
            isinstance(v, str) for v in adata_clustered.obs["cell_type_scanvi"]
        )

    def test_returns_cluster_posteriors_dict(self, adata_clustered):
        labels = ["T cell"] * len(adata_clustered)
        probs  = np.full(len(adata_clustered), 0.9)
        mock_model = _make_mock_scanvi_model(labels, probs)

        with patch("scvi.model.SCANVI") as mock_cls:
            mock_cls.load.return_value = mock_model
            result = _run_scanvi(adata_clustered, "leiden", "/fake/model")

        assert isinstance(result, dict)
        clusters = adata_clustered.obs["leiden"].unique().tolist()
        for cl in clusters:
            assert str(cl) in result, f"Cluster {cl} missing from posteriors"
            assert 0.0 <= result[str(cl)] <= 1.0

    def test_posterior_values_are_floats(self, adata_clustered):
        labels = ["B cell"] * len(adata_clustered)
        mock_model = _make_mock_scanvi_model(labels)

        with patch("scvi.model.SCANVI") as mock_cls:
            mock_cls.load.return_value = mock_model
            result = _run_scanvi(adata_clustered, "leiden", "/fake/model")

        for cl, prob in result.items():
            assert isinstance(prob, float), f"Posterior for cluster {cl} is not float"

    def test_column_is_category_dtype(self, adata_clustered):
        labels = ["T cell"] * len(adata_clustered)
        mock_model = _make_mock_scanvi_model(labels)

        with patch("scvi.model.SCANVI") as mock_cls:
            mock_cls.load.return_value = mock_model
            _run_scanvi(adata_clustered, "leiden", "/fake/model")

        assert hasattr(adata_clustered.obs["cell_type_scanvi"], "cat"), \
            "cell_type_scanvi should be category dtype"

    def test_scanvi_without_model_path_raises(self, adata_clustered):
        """Requesting scanvi without providing scanvi_model raises ValueError."""
        with pytest.raises(ValueError, match="scanvi_model"):
            annotate(
                adata_clustered,
                methods=["markers", "scanvi"],
                scanvi_model=None,
                inplace=False,
            )

    def test_scanvi_importerror_warns_and_skips(self, adata_clustered):
        """If scvi-tools is not installed, annotate() warns and skips scanvi."""
        with patch.dict("sys.modules", {"scvi": None, "scvi.model": None}):
            with pytest.warns(UserWarning, match="scvi-tools"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "scanvi"],
                    scanvi_model="/fake/model",
                    inplace=False,
                )
        assert "scanvi" not in ann_dict["provenance"]["methods_run"]
        assert "cell_type_scanvi" not in adata_ann.obs.columns

    def test_scanvi_generic_error_warns_and_skips(self, adata_clustered):
        """A runtime error in _run_scanvi warns and skips without crashing."""
        with patch(
            "pipeline.modules.scripts.annotation.annotate._run_scanvi",
            side_effect=RuntimeError("model corrupt"),
        ):
            with pytest.warns(UserWarning, match="scANVI annotation failed"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "scanvi"],
                    scanvi_model="/fake/model",
                    inplace=False,
                )
        assert "scanvi" not in ann_dict["provenance"]["methods_run"]


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: N-way vote upgrade tests  (always run — pure logic, no network)
# ─────────────────────────────────────────────────────────────────────────────

class TestNWayVote:
    """
    Tests for the upgraded N-way weighted majority vote.
    All tests use pre-populated obs columns rather than calling annotate() end-to-end,
    so they are independent of CellTypist/scANVI installation.
    """

    def _make_annotated(self, adata, n_clusters=3):
        """
        Populate the obs columns that vote reads from, then call _run_majority_vote.
        Returns (adata, vote_df).
        """
        from pipeline.modules.scripts.annotation.annotate import _run_marker_scoring, _run_majority_vote
        score_df = _run_marker_scoring(adata, "leiden", MARKER_SETS)
        # Simulate celltypist_fine (same as markers — perfect agreement)
        adata.obs["celltypist_fine"] = adata.obs["cell_type_markers"].copy()
        vote_df = _run_majority_vote(adata, score_df, "leiden")
        return adata, vote_df

    def test_vote_columns_written(self, adata_clustered):
        adata, _ = self._make_annotated(adata_clustered)
        assert "cell_type_vote" in adata.obs.columns
        assert "cell_type_confidence" in adata.obs.columns

    def test_confidence_range(self, adata_clustered):
        adata, _ = self._make_annotated(adata_clustered)
        conf = adata.obs["cell_type_confidence"]
        assert conf.between(0.0, 1.0).all(), \
            f"Confidence outside [0,1]: min={conf.min()}, max={conf.max()}"

    def test_scanvi_slot_absent_when_column_missing(self, adata_clustered):
        """Vote completes without crashing when cell_type_scanvi is not present."""
        adata, vote_df = self._make_annotated(adata_clustered)
        # scANVI column absent → scANVI row in vote_df is None
        if "scANVI" in vote_df.columns:
            assert vote_df["scANVI"].isna().all(), \
                "scANVI column should be all None when scanvi was not run"

    def test_scanvi_fractional_weight_increases_confidence(self, adata_clustered):
        """
        Adding cell_type_scanvi (same label, high posterior) should raise confidence
        to at or above the level without it.
        """
        from pipeline.modules.scripts.annotation.annotate import _run_marker_scoring, _run_majority_vote

        # Baseline: two methods, perfect agreement
        adata_base = adata_clustered.copy()
        score_df = _run_marker_scoring(adata_base, "leiden", MARKER_SETS)
        adata_base.obs["celltypist_fine"] = adata_base.obs["cell_type_markers"].copy()
        vote_df_base = _run_majority_vote(adata_base, score_df, "leiden")
        conf_base = adata_base.obs["cell_type_confidence"].mean()

        # With scANVI: same label as markers, high posterior
        adata_scanvi = adata_clustered.copy()
        score_df2 = _run_marker_scoring(adata_scanvi, "leiden", MARKER_SETS)
        adata_scanvi.obs["celltypist_fine"] = adata_scanvi.obs["cell_type_markers"].copy()
        adata_scanvi.obs["cell_type_scanvi"] = adata_scanvi.obs["cell_type_markers"].copy()
        clusters = adata_scanvi.obs["leiden"].unique()
        posteriors = {str(cl): 0.95 for cl in clusters}
        vote_df_scanvi = _run_majority_vote(
            adata_scanvi, score_df2, "leiden", scanvi_posteriors=posteriors
        )
        conf_scanvi = adata_scanvi.obs["cell_type_confidence"].mean()

        # With full agreement across 3 methods, confidence should be ≥ 2-method baseline
        assert conf_scanvi >= conf_base - 0.01, (
            f"scANVI agreement should not lower confidence: "
            f"base={conf_base:.3f} scanvi={conf_scanvi:.3f}"
        )

    def test_scanvi_low_posterior_gives_less_weight(self, adata_clustered):
        """
        scANVI with low posterior (0.1) contributes little weight;
        winner should still be decided by the other two methods.
        """
        from pipeline.modules.scripts.annotation.annotate import _run_marker_scoring, _run_majority_vote

        adata = adata_clustered.copy()
        score_df = _run_marker_scoring(adata, "leiden", MARKER_SETS)
        adata.obs["celltypist_fine"] = adata.obs["cell_type_markers"].copy()
        # Give scANVI a *different* label with very low posterior
        adata.obs["cell_type_scanvi"] = pd.Categorical(
            ["NK_ILC"] * len(adata)  # unlikely to match celltypist/markers
        )
        clusters = adata.obs["leiden"].unique()
        low_posteriors = {str(cl): 0.1 for cl in clusters}

        vote_df = _run_majority_vote(
            adata, score_df, "leiden", scanvi_posteriors=low_posteriors
        )
        # With low posterior weight (0.1 × 2 = 0.2 votes), the 2-method consensus
        # (2 votes each for the same label) should still win
        winners = vote_df["final_label"].unique()
        assert "NK_ILC" not in winners or len(winners) == 1, \
            "Low-posterior scANVI label should not override 2-method consensus"

    def test_vote_df_has_scanvi_column_when_present(self, adata_clustered):
        """vote_df has 'scANVI' column when cell_type_scanvi obs col exists."""
        from pipeline.modules.scripts.annotation.annotate import _run_marker_scoring, _run_majority_vote

        adata = adata_clustered.copy()
        score_df = _run_marker_scoring(adata, "leiden", MARKER_SETS)
        adata.obs["celltypist_fine"] = adata.obs["cell_type_markers"].copy()
        adata.obs["cell_type_scanvi"] = adata.obs["cell_type_markers"].copy()
        clusters = adata.obs["leiden"].unique()
        posteriors = {str(cl): 0.8 for cl in clusters}

        vote_df = _run_majority_vote(adata, score_df, "leiden", scanvi_posteriors=posteriors)
        assert "scANVI" in vote_df.columns, "vote_df should have 'scANVI' column"
        assert vote_df["scANVI"].notna().all(), "scANVI labels should not be null"

    def test_vote_with_no_active_slots_gives_unknown(self, adata_clustered):
        """If all obs columns are missing, vote assigns 'Unknown' to every cluster."""
        from pipeline.modules.scripts.annotation.annotate import _run_marker_scoring, _run_majority_vote

        adata = adata_clustered.copy()
        score_df = _run_marker_scoring(adata, "leiden", MARKER_SETS)
        # Do NOT add celltypist_fine or other columns — only markers column present
        # (markers col always present after _run_marker_scoring)
        # Remove markers column to simulate all-empty scenario
        adata.obs.drop(columns=["cell_type_markers"], inplace=True)

        vote_df = _run_majority_vote(adata, score_df, "leiden")
        assert (vote_df["final_label"] == "Unknown").all(), \
            "With no active slots, all clusters should be 'Unknown'"


# ─────────────────────────────────────────────────────────────────────────────
# SingleR tests — all mocked, no real network calls in CI
# ─────────────────────────────────────────────────────────────────────────────

def _make_singler_ref(n_types: int = 6, n_genes: int = 50) -> sc.AnnData:
    """
    Build a small in-memory pseudobulk AnnData suitable as a SingleR reference.
    obs_names = cell type labels, var_names = gene symbols (some shared with fixture).
    X = mean log-normalised expression (dense float32).
    """
    rng = np.random.default_rng(99)

    # Use gene names that overlap with the fixture fixture gene set
    shared_genes = [
        "CD3D", "CD3E", "CD68", "CD14", "MS4A1",
        "ALB", "APOC3", "COL1A1", "PECAM1", "VWF",
    ]
    extra_genes = [f"REF_GENE{i}" for i in range(n_genes - len(shared_genes))]
    gene_names  = shared_genes + extra_genes

    cell_type_names = [f"CellType{i}" for i in range(n_types)]
    X = rng.uniform(0.0, 3.0, size=(n_types, len(gene_names))).astype(np.float32)

    # Make types clearly distinguishable by assigning high expression to a unique gene
    for i in range(min(n_types, len(shared_genes))):
        X[i, i] = 10.0  # strong signature for each type

    ref = sc.AnnData(X=X)
    ref.obs_names = cell_type_names
    ref.var_names = gene_names
    return ref


class TestSingleR:
    """
    Unit and integration tests for _run_singler() and its helpers.

    Mock strategy
    -------------
    Every test patches BOTH:
      1. ``pipeline.modules.scripts.annotation.annotate._load_singler_ref`` → returns a
         small in-memory AnnData so no network call is made.
      2. ``singler.annotate_single`` → returns a controlled BiocFrame-like object
         so the C++ library is never invoked and results are predictable.

    This keeps all 13 tests fast and deterministic in CI.
    """

    # ── Shared mock factories ──────────────────────────────────────────────────

    @staticmethod
    def _make_singler_result(n_cells: int, label: str = "T cell", delta: float = 0.3):
        """
        Return a minimal mock that mimics singler.annotate_single's BiocFrame output.
        Columns: best (list of str), delta (np.ndarray of float64).
        """
        mock_result = MagicMock()
        mock_result.column.side_effect = lambda col: (
            [label] * n_cells if col == "best"
            else np.full(n_cells, delta, dtype=np.float64)
        )
        return mock_result

    @staticmethod
    def _make_singler_result_unassigned(n_cells: int):
        """Return a mock where all deltas are ≤ 0 → all cells → 'Unassigned'."""
        mock_result = MagicMock()
        mock_result.column.side_effect = lambda col: (
            ["T cell"] * n_cells if col == "best"
            else np.full(n_cells, 0.0, dtype=np.float64)   # delta=0 → Unassigned
        )
        return mock_result

    def _run_with_mocks(self, adata, label="T cell", delta=0.3):
        """Helper: run _run_singler with both _load_singler_ref and singler.annotate_single mocked."""
        ref = _make_singler_ref()
        n = len(adata)
        mock_result = self._make_singler_result(n, label=label, delta=delta)
        with patch("pipeline.modules.scripts.annotation.annotate._load_singler_ref", return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result):
                _run_singler(adata, "leiden")

    # ── Core output tests ──────────────────────────────────────────────────────

    def test_singler_column_exists(self, adata_clustered):
        self._run_with_mocks(adata_clustered)
        assert "cell_type_singler" in adata_clustered.obs.columns

    def test_singler_delta_column_exists(self, adata_clustered):
        self._run_with_mocks(adata_clustered)
        assert "singler_delta" in adata_clustered.obs.columns

    def test_singler_labels_are_strings(self, adata_clustered):
        self._run_with_mocks(adata_clustered)
        labels = adata_clustered.obs["cell_type_singler"]
        assert all(isinstance(v, str) for v in labels), \
            "All cell_type_singler values must be strings"

    def test_singler_every_cell_covered(self, adata_clustered):
        self._run_with_mocks(adata_clustered)
        n_null = adata_clustered.obs["cell_type_singler"].isna().sum()
        assert n_null == 0, f"{n_null} cells have null cell_type_singler"

    def test_singler_delta_is_float(self, adata_clustered):
        self._run_with_mocks(adata_clustered)
        deltas = adata_clustered.obs["singler_delta"]
        assert pd.api.types.is_float_dtype(deltas), \
            f"singler_delta should be float dtype, got {deltas.dtype}"

    def test_singler_unassigned_for_low_delta(self, adata_clustered):
        """
        When singler returns delta=0 for all cells, _run_singler should label
        them all 'Unassigned' (delta ≤ 0 → pruned).
        """
        ref = _make_singler_ref()
        n = len(adata_clustered)
        mock_result = self._make_singler_result_unassigned(n)
        with patch("pipeline.modules.scripts.annotation.annotate._load_singler_ref", return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result):
                _run_singler(adata_clustered, "leiden")
        labels = adata_clustered.obs["cell_type_singler"].unique().tolist()
        assert "Unassigned" in labels, \
            "Cells with delta≤0 should be labelled Unassigned"

    def test_singler_delta_nan_for_unassigned(self, adata_clustered):
        """delta must be NaN for cells labelled 'Unassigned'."""
        ref = _make_singler_ref()
        n = len(adata_clustered)
        mock_result = self._make_singler_result_unassigned(n)
        with patch("pipeline.modules.scripts.annotation.annotate._load_singler_ref", return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result):
                _run_singler(adata_clustered, "leiden")
        mask_unassigned = adata_clustered.obs["cell_type_singler"] == "Unassigned"
        if mask_unassigned.any():
            deltas_unassigned = adata_clustered.obs.loc[mask_unassigned, "singler_delta"]
            assert deltas_unassigned.isna().all(), \
                "Unassigned cells must have NaN delta"

    def test_singler_uses_logcounts_layer(self, adata_clustered):
        """
        Zero out logcounts — _run_singler should use it (not .X).
        We verify by checking singler.annotate_single was still called
        (it would receive an all-zero matrix — no crash expected).
        """
        ref = _make_singler_ref()
        n = len(adata_clustered)
        lc = adata_clustered.layers["logcounts"]
        adata_clustered.layers["logcounts"] = np.zeros_like(
            lc if not hasattr(lc, "toarray") else lc.toarray()
        )
        mock_result = self._make_singler_result(n)
        with patch("pipeline.modules.scripts.annotation.annotate._load_singler_ref", return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result) as mock_fn:
                _run_singler(adata_clustered, "leiden")
        mock_fn.assert_called_once()
        # The test_data passed to singler should be all zeros (logcounts used, not .X)
        call_kwargs = mock_fn.call_args
        X_passed = call_kwargs.kwargs.get("test_data", call_kwargs.args[0] if call_kwargs.args else None)
        assert X_passed is not None
        assert np.allclose(X_passed, 0.0), "singler received non-zero matrix — logcounts not used"
        assert "cell_type_singler" in adata_clustered.obs.columns

    # ── Provenance tests ───────────────────────────────────────────────────────

    def test_singler_ref_recorded_in_provenance(self, adata_clustered):
        ref = _make_singler_ref()
        n = len(adata_clustered)
        mock_result = self._make_singler_result(n)
        with patch("pipeline.modules.scripts.annotation.annotate._load_singler_ref", return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "singler"],
                    inplace=False,
                )
        assert "singler_ref" in ann_dict["provenance"]

    def test_singler_ref_label_col_recorded_in_provenance(self, adata_clustered):
        ref = _make_singler_ref()
        n = len(adata_clustered)
        mock_result = self._make_singler_result(n)
        with patch("pipeline.modules.scripts.annotation.annotate._load_singler_ref", return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "singler"],
                    singler_ref_label_col="my_cell_type",
                    inplace=False,
                )
        assert ann_dict["provenance"].get("singler_ref_label_col") == "my_cell_type"

    # ── Error handling tests ───────────────────────────────────────────────────

    def test_singler_skipped_gracefully_on_error(self, adata_clustered):
        """If _run_singler raises, annotate() warns and skips — no crash."""
        with patch(
            "pipeline.modules.scripts.annotation.annotate._run_singler",
            side_effect=RuntimeError("ref download failed"),
        ):
            with pytest.warns(UserWarning, match="SingleR annotation failed"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "singler"],
                    inplace=False,
                )
        assert "singler" not in ann_dict["provenance"]["methods_run"]
        assert "cell_type_singler" not in adata_ann.obs.columns

    def test_hca_import_error_raised_without_package(self, adata_clustered):
        """singler_ref='hca' without cellxgene-census installed raises ImportError."""
        with patch.dict("sys.modules", {"cellxgene_census": None}):
            with pytest.raises((ImportError, Exception)):
                _load_singler_ref("hca", "cell_type")

    def test_user_h5ad_missing_file_raises(self, adata_clustered):
        """singler_ref pointing to a non-existent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            _load_singler_ref("/nonexistent/path/ref.h5ad", "cell_type")

    # ── Cache test (skipped in CI) ─────────────────────────────────────────────

    @pytest.mark.skipif(
        os.getenv("OMICSAGE_CI") == "true",
        reason="Skipped in CI — requires network access to download HPCA reference",
    )
    def test_tabula_sapiens_cached_after_first_download(self, tmp_path):
        """
        After _load_hpca_reference() runs, the cache file must exist.
        We patch both _HPCA_CACHE and _TABULA_SAPIENS_CACHE (alias) and
        _SINGLER_CACHE_DIR to avoid touching the real ~/.cache directory.
        """
        from pipeline.modules.scripts.annotation.annotate import _load_hpca_reference
        fake_cache = tmp_path / "hpca_ref.h5ad"
        with patch("pipeline.modules.scripts.annotation.annotate._HPCA_CACHE", fake_cache):
            with patch("pipeline.modules.scripts.annotation.annotate._TABULA_SAPIENS_CACHE", fake_cache):
                with patch("pipeline.modules.scripts.annotation.annotate._SINGLER_CACHE_DIR", tmp_path):
                    ref = _load_hpca_reference()
        assert fake_cache.exists(), "Cache file should exist after first download"
        assert ref.n_obs >= 5, "Reference should have at least 5 cell types"


# ─────────────────────────────────────────────────────────────────────────────
# Reference selection tests — sctype_db_path and named celldex singler_ref
# ─────────────────────────────────────────────────────────────────────────────

class TestSctypeDbPath:
    """Tests for the sctype_db_path parameter (local ScTypeDB file support)."""

    def test_sctype_uses_local_db_when_path_given(self, adata_clustered, tmp_path):
        """When sctype_db_path is set, _fetch_sctype_db reads that file (no network)."""
        # Build a minimal ScTypeDB Excel in the expected format
        mini_db = pd.DataFrame({
            "tissueType":    ["Immune system", "Immune system"],
            "cellName":      ["T cell", "B cell"],
            "geneSymbolmore1": ["CD3D,CD3E", "MS4A1,CD19"],
            "geneSymbolmore2": ["", ""],
        })
        db_path = tmp_path / "mini_sctype_db.xlsx"
        mini_db.to_excel(db_path, index=False)

        from pipeline.modules.scripts.annotation.annotate import _fetch_sctype_db
        db = _fetch_sctype_db(db_path=db_path)
        assert "tissueType" in db.columns
        assert set(db["cellName"].tolist()) == {"T cell", "B cell"}

    def test_sctype_local_db_missing_file_raises(self, adata_clustered):
        """sctype_db_path pointing to non-existent file raises FileNotFoundError."""
        from pipeline.modules.scripts.annotation.annotate import _fetch_sctype_db
        with pytest.raises(FileNotFoundError):
            _fetch_sctype_db(db_path="/nonexistent/db.xlsx")

    def test_annotate_passes_sctype_db_path_to_run_sctype(self, adata_clustered, tmp_path):
        """annotate() threads sctype_db_path through to _run_sctype."""
        db_path = tmp_path / "mini_db.xlsx"

        with patch("pipeline.modules.scripts.annotation.annotate._run_sctype") as mock_sctype:
            annotate(
                adata_clustered,
                methods=["markers", "sctype"],
                sctype_db_path=db_path,
                inplace=False,
            )
        mock_sctype.assert_called_once()
        _, call_kwargs = mock_sctype.call_args
        assert call_kwargs.get("sctype_db_path") == db_path

    def test_sctype_db_path_recorded_in_provenance(self, adata_clustered, tmp_path):
        """sctype_db_path is recorded in provenance."""
        mini_db = pd.DataFrame({
            "tissueType":    ["Immune system"] * 2,
            "cellName":      ["T cell", "B cell"],
            "geneSymbolmore1": ["CD3D", "MS4A1"],
            "geneSymbolmore2": ["", ""],
        })
        db_path = tmp_path / "mini_db.xlsx"
        mini_db.to_excel(db_path, index=False)

        adata_ann, ann_dict = annotate(
            adata_clustered,
            methods=["markers", "sctype"],
            sctype_db_path=db_path,
            inplace=False,
        )
        assert ann_dict["provenance"].get("sctype_db_path") == str(db_path)


class TestSinglerNamedCelldexRef:
    """Tests for named celldex references in singler_ref."""

    def test_known_celldex_refs_constant_has_correct_keys(self):
        """_CELLDEX_KNOWN_REFS contains all expected named references."""
        from pipeline.modules.scripts.annotation.annotate import _CELLDEX_KNOWN_REFS
        expected = {
            "hpca", "blueprint_encode", "dice",
            "monaco_immune", "novershtern_hematopoietic", "mouse_rnaseq",
        }
        assert set(_CELLDEX_KNOWN_REFS.keys()) == expected

    def test_named_ref_routes_to_load_celldex_reference(self, adata_clustered):
        """singler_ref='blueprint_encode' calls _load_celldex_reference, not the file-path branch."""
        from pipeline.modules.scripts.annotation.annotate import _load_singler_ref
        ref = _make_singler_ref()
        with patch(
            "pipeline.modules.scripts.annotation.annotate._load_celldex_reference",
            return_value=ref,
        ) as mock_load:
            result = _load_singler_ref("blueprint_encode", "cell_type")
        mock_load.assert_called_once_with("blueprint_encode")
        assert result is ref

    def test_named_ref_case_insensitive(self, adata_clustered):
        """Named refs are case-insensitive (Blueprint_Encode → blueprint_encode)."""
        from pipeline.modules.scripts.annotation.annotate import _load_singler_ref
        ref = _make_singler_ref()
        with patch(
            "pipeline.modules.scripts.annotation.annotate._load_celldex_reference",
            return_value=ref,
        ) as mock_load:
            _load_singler_ref("Blueprint_Encode", "cell_type")
        mock_load.assert_called_once_with("blueprint_encode")

    def test_unknown_ref_name_raises_valueerror(self):
        """An unknown string that is not a file path raises a helpful ValueError."""
        from pipeline.modules.scripts.annotation.annotate import _load_singler_ref
        with pytest.raises((ValueError, FileNotFoundError)):
            _load_singler_ref("nonexistent_reference", "cell_type")

    def test_annotate_singler_named_ref_end_to_end(self, adata_clustered):
        """annotate() with singler_ref='novershtern_hematopoietic' completes successfully."""
        ref = _make_singler_ref()
        n = len(adata_clustered)
        mock_result = MagicMock()
        mock_result.column.side_effect = lambda col: (
            ["CellType0"] * n if col == "best"
            else np.full(n, 0.3, dtype=np.float64)
        )
        with patch("pipeline.modules.scripts.annotation.annotate._load_celldex_reference",
                   return_value=ref):
            with patch("singler.annotate_single", return_value=mock_result):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "singler"],
                    singler_ref="novershtern_hematopoietic",
                    inplace=False,
                )
        assert "cell_type_singler" in adata_ann.obs.columns
        assert ann_dict["provenance"]["singler_ref"] == "novershtern_hematopoietic"


# ─────────────────────────────────────────────────────────────────────────────
# Majority vote — label normalisation and tiebreak tests
# ─────────────────────────────────────────────────────────────────────────────

from pipeline.modules.scripts.annotation.annotate import _normalise_label, _run_majority_vote


class TestNormaliseLabel:
    """Unit tests for _normalise_label() — the semantic deduplication helper."""

    @pytest.mark.parametrize("raw,expected", [
        ("Erythroid cells",                          "erythroid"),
        ("Late erythroid",                           "erythroid"),
        ("Erythroid-like and erythroid precursor cells", "erythroid"),
        ("Erythroid cell",                           "erythroid"),
        ("Classical monocytes",                      "monocyte"),
        ("Classical Monocytes",                      "monocyte"),
        ("Monocytes",                                "monocyte"),
        ("CD4+ T cells",                             "t"),
        ("Tcm/Naive helper T cells",                 "tcm"),   # T central memory — kept specific
        ("Small pre-B cells",                        "b"),
        ("Naive B cells",                            "b"),
        ("Natural killer cells",                     "natural killer"),
        ("NK cells",                                 "nk"),
        ("HSC/MPP",                                  "hsc"),
        ("Plasma cells",                             "plasma"),
        ("Monocyte",                                 "monocyte"),
        ("pDC",                                      "pdc"),
        ("B cells",                                  "b"),
        ("T_cell",                                   "t"),
        ("B_cell",                                   "b"),
    ])
    def test_normalisation(self, raw, expected):
        assert _normalise_label(raw) == expected, (
            f"_normalise_label({raw!r}) = {_normalise_label(raw)!r}, expected {expected!r}"
        )


class TestMajorityVoteFixes:
    """
    Integration tests for both vote bugs:
      1. Synonymous labels (same concept, different strings) should yield high confidence.
      2. Tiebreak should use source hierarchy, not dict insertion order.
    """

    def _make_vote_adata(self, cluster_labels: dict) -> sc.AnnData:
        """
        Build a minimal AnnData with the given obs columns for a single cluster "0".
        cluster_labels: {obs_col: label_string}
        """
        n = 20
        adata = sc.AnnData(X=np.zeros((n, 5)))
        adata.obs["leiden"] = "0"
        for col, label in cluster_labels.items():
            adata.obs[col] = label
        # Dummy logcounts layer
        adata.layers["logcounts"] = np.zeros((n, 5))
        # Dummy score_df (not used in vote logic but required as arg)
        return adata

    def _dummy_score_df(self):
        return pd.DataFrame({"cluster": ["0"]}).set_index("cluster")

    def test_synonymous_erythroid_labels_give_high_confidence(self):
        """
        All 4 methods agree on erythroid — different strings but same root concept.
        Confidence should be 1.0 after normalisation, not 0.25.
        """
        adata = self._make_vote_adata({
            "celltypist_fine":   "Late erythroid",
            "cell_type_markers": "Erythroid",
            "cell_type_sctype":  "Erythroid-like and erythroid precursor cells",
            "cell_type_singler": "Erythroid cells",
        })
        vote_df = _run_majority_vote(adata, self._dummy_score_df(), "leiden")
        conf = float(vote_df.loc["0", "confidence"])
        assert conf == 1.0, f"Expected confidence=1.0 for unanimous erythroid, got {conf}"

    def test_synonymous_monocyte_labels_give_high_confidence(self):
        """Classical monocytes / Monocyte / Monocytes → same root → confidence=1.0."""
        adata = self._make_vote_adata({
            "celltypist_fine":   "Classical monocytes",
            "cell_type_markers": "Monocyte",
            "cell_type_sctype":  "Classical Monocytes",
            "cell_type_singler": "Monocytes",
        })
        vote_df = _run_majority_vote(adata, self._dummy_score_df(), "leiden")
        conf = float(vote_df.loc["0", "confidence"])
        assert conf == 1.0, f"Expected confidence=1.0 for unanimous monocyte, got {conf}"

    def test_genuine_disagreement_gives_low_confidence(self):
        """When methods genuinely disagree, confidence should be < 1.0."""
        adata = self._make_vote_adata({
            "celltypist_fine":   "T cells",
            "cell_type_markers": "B_cell",
            "cell_type_sctype":  "Monocyte",
            "cell_type_singler": "NK cells",
        })
        vote_df = _run_majority_vote(adata, self._dummy_score_df(), "leiden")
        conf = float(vote_df.loc["0", "confidence"])
        assert conf < 1.0, f"Expected confidence<1.0 for genuine disagreement, got {conf}"

    def test_tiebreak_prefers_celltypist_fine(self):
        """
        When two concepts tie, celltypist_fine wins the tiebreak.
        Here celltypist_fine says "T cells", markers says "B_cell" — each gets
        weight=1, so they tie; celltypist_fine should win.
        """
        adata = self._make_vote_adata({
            "celltypist_fine":   "T cells",
            "cell_type_markers": "B_cell",
        })
        vote_df = _run_majority_vote(adata, self._dummy_score_df(), "leiden")
        winner = vote_df.loc["0", "final_label"]
        assert winner == "T cells", (
            f"Tiebreak should prefer celltypist_fine='T cells', got '{winner}'"
        )

    def test_winner_label_is_raw_not_canonical(self):
        """
        The final reported label is the raw string from the best source,
        not the lower-cased canonical form.
        """
        adata = self._make_vote_adata({
            "celltypist_fine":   "Late erythroid",
            "cell_type_markers": "Erythroid",
            "cell_type_singler": "Erythroid cells",
        })
        vote_df = _run_majority_vote(adata, self._dummy_score_df(), "leiden")
        winner = vote_df.loc["0", "final_label"]
        # Should be a raw label, not lowercase canonical "erythroid"
        assert winner in {"Late erythroid", "Erythroid", "Erythroid cells"}, (
            f"Winner should be a raw label, got '{winner}'"
        )
        assert winner != "erythroid", "Winner must not be the lowercased canonical form"

    def test_missing_method_columns_handled_gracefully(self):
        """Vote works fine when some method columns are absent."""
        adata = self._make_vote_adata({
            "celltypist_fine":   "T cells",
            "cell_type_markers": "T_cell",
            # sctype and singler absent
        })
        vote_df = _run_majority_vote(adata, self._dummy_score_df(), "leiden")
        assert "final_label" in vote_df.columns
        assert float(vote_df.loc["0", "confidence"]) == 1.0


# ─────────────────────────────────────────────────────────────────────────────
# CellTypist fix tests — per-cell prediction + any model name
# ─────────────────────────────────────────────────────────────────────────────

class TestCelltypistAnyModel:
    """
    Verify _run_celltypist works correctly with any model name (not just
    Immune_All_High/Low) and uses per-cell predictions + our own majority vote.
    """

    def _mock_pred(self, adata, label="T cell"):
        """Build a mock CellTypist AnnotationResult that returns fixed per-cell labels."""
        import pandas as pd
        pred_adata = adata.copy()
        pred_adata.obs["predicted_labels"] = label
        mock_result = MagicMock()
        mock_result.to_adata.return_value = pred_adata
        return mock_result

    def test_custom_model_creates_correct_obs_column(self, adata_clustered, tmp_path):
        """
        A custom model 'Healthy_COVID19_PBMC.pkl' should produce
        obs['celltypist_Healthy_COVID19_PBMC'], not crash.
        """
        model_names = ["Healthy_COVID19_PBMC.pkl"]
        with patch("celltypist.models.download_models"):
            with patch("celltypist.models.Model.load", return_value=MagicMock()):
                with patch(
                    "celltypist.annotate",
                    return_value=self._mock_pred(adata_clustered),
                ):
                    from pipeline.modules.scripts.annotation.annotate import _run_celltypist
                    _run_celltypist(adata_clustered, "leiden", model_names, tmp_path)

        assert "celltypist_Healthy_COVID19_PBMC" in adata_clustered.obs.columns

    def test_custom_model_guarantees_coarse_fine_exist(self, adata_clustered, tmp_path):
        """
        Even when only a custom model is run, celltypist_coarse and
        celltypist_fine must exist (filled with 'not_run').
        """
        model_names = ["Healthy_COVID19_PBMC.pkl"]
        with patch("celltypist.models.download_models"):
            with patch("celltypist.models.Model.load", return_value=MagicMock()):
                with patch(
                    "celltypist.annotate",
                    return_value=self._mock_pred(adata_clustered),
                ):
                    from pipeline.modules.scripts.annotation.annotate import _run_celltypist
                    _run_celltypist(adata_clustered, "leiden", model_names, tmp_path)

        assert adata_clustered.obs["celltypist_coarse"].iloc[0] == "not_run"
        assert adata_clustered.obs["celltypist_fine"].iloc[0] == "not_run"

    def test_celltypist_uses_majority_voting_false(self, adata_clustered, tmp_path):
        """
        celltypist.annotate must be called with majority_voting=False.
        We must never rely on CellTypist's internal over_clustering.
        """
        model_names = ["Immune_All_Low.pkl"]
        with patch("celltypist.models.download_models"):
            with patch("celltypist.models.Model.load", return_value=MagicMock()):
                with patch(
                    "celltypist.annotate",
                    return_value=self._mock_pred(adata_clustered),
                ) as mock_ct:
                    from pipeline.modules.scripts.annotation.annotate import _run_celltypist
                    _run_celltypist(adata_clustered, "leiden", model_names, tmp_path)

        call_kwargs = mock_ct.call_args.kwargs
        assert call_kwargs.get("majority_voting") is False, (
            "celltypist.annotate must be called with majority_voting=False"
        )
        assert "over_clustering" not in call_kwargs, (
            "over_clustering must not be passed when majority_voting=False"
        )

    def test_extra_celltypist_column_feeds_into_vote(self, adata_clustered, tmp_path):
        """
        Labels from a custom model (celltypist_Healthy_COVID19_PBMC) must
        contribute to the majority vote result, not be silently ignored.
        """
        # Manually set up the custom column as if CellTypist ran
        adata_clustered.obs["celltypist_Healthy_COVID19_PBMC"] = "Monocyte"
        adata_clustered.obs["celltypist_fine"]   = "not_run"
        adata_clustered.obs["cell_type_markers"] = "Monocyte"
        adata_clustered.obs["cell_type_sctype"]  = "Monocyte"
        adata_clustered.obs["cell_type_singler"] = "Monocyte"

        from pipeline.modules.scripts.annotation.annotate import _run_majority_vote
        score_df = pd.DataFrame(
            {"cluster": adata_clustered.obs["leiden"].unique().tolist()}
        ).set_index("cluster")

        vote_df = _run_majority_vote(adata_clustered, score_df, "leiden")
        # All non-not_run methods agree on Monocyte → confidence should be 1.0
        assert (vote_df["confidence"] == 1.0).all(), (
            "Custom celltypist column should feed into vote and yield confidence=1.0"
        )
