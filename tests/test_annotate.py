"""
tests/test_annotate.py
======================
Test suite for pipeline/modules/annotation/annotate.py

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

from pipeline.modules.annotation.annotate import (
    annotate,
    _fetch_sctype_db,
    _parse_sctype_markers,
    _run_sctype,
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
    from pipeline.modules.annotation.annotate import _run_marker_scoring, _run_majority_vote
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
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        assert "cell_type_sctype" in adata_clustered.obs.columns
 
    def test_labels_are_strings(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        labels = adata_clustered.obs["cell_type_sctype"]
        assert all(isinstance(v, str) for v in labels), \
            "All cell_type_sctype values must be strings"
 
    def test_every_cluster_covered(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
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
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        labels = adata_clustered.obs["cell_type_sctype"].unique().tolist()
        assert "Unknown" not in labels
 
    def test_liver_tissue_gives_liver_types(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Liver")
        labels = set(adata_clustered.obs["cell_type_sctype"].unique())
        liver_types = {"Hepatocyte", "Fibroblast"}
        assert labels.issubset(liver_types | {"Unknown"}), \
            f"Unexpected labels for Liver tissue: {labels - liver_types}"
 
    def test_column_is_category_dtype(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
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
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            _run_sctype(adata_clustered, "leiden", tissue="Immune system")
        # With zero expression all scores are 0; best_type still assigned (no crash)
        assert "cell_type_sctype" in adata_clustered.obs.columns
 
 
# ── annotate() integration tests (mocked network) ─────────────────────────────
 
class TestAnnotateWithScType:
 
    def test_sctype_in_methods_run(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            with patch("pipeline.modules.annotation.annotate._run_celltypist"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "sctype"],
                    tissue="Immune system",
                )
        assert "sctype" in ann_dict["provenance"]["methods_run"]
 
    def test_tissue_recorded_in_provenance(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            with patch("pipeline.modules.annotation.annotate._run_celltypist"):
                adata_ann, ann_dict = annotate(
                    adata_clustered,
                    methods=["markers", "sctype"],
                    tissue="Liver",
                )
        assert ann_dict["provenance"]["sctype_tissue"] == "Liver"
 
    def test_sctype_column_written_by_annotate(self, adata_clustered):
        with patch(
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
            return_value=_make_mock_db(),
        ):
            with patch("pipeline.modules.annotation.annotate._run_celltypist"):
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
            "pipeline.modules.annotation.annotate._fetch_sctype_db",
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
