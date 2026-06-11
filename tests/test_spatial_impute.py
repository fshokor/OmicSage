"""
tests/test_spatial_impute.py — OmicSage Phase 7 extension

Tests for pipeline/modules/spatial/spatial_impute.py and
reports/templates/spatial/spatial_impute_report.py.

Strategy
--------
- All tests use tiny synthetic data (20 spots × 50 genes, 100 cells × 50 genes)
  so no real datasets or GPU are required.
- Tangram is mocked throughout — `pytest.importorskip` pattern is NOT used
  here because we need to test even when tangram-sc is NOT installed (CI).
  Instead, we patch the module-level flag `_TANGRAM_AVAILABLE` and inject a
  mock `tg` namespace for every test that exercises the Tangram code path.
- gimVI path is tested only for its ImportError guard (scvi-tools not required).
- Report rendering is tested on mock data; no squidpy spatial plots are
  attempted (squidpy is skipped gracefully inside the report).
"""

from __future__ import annotations

import importlib
import sys
import types
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def rng():
    return np.random.default_rng(42)


@pytest.fixture()
def spatial_adata(rng) -> ad.AnnData:
    """20 spots × 50 genes with raw counts + spatial coordinates."""
    n_spots, n_genes = 20, 50
    X = rng.integers(0, 50, size=(n_spots, n_genes)).astype(np.float32)
    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
    obs_names  = [f"spot_{i}" for i in range(n_spots)]
    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=gene_names),
    )
    # Spatial coordinates (required by sq.pl.spatial_scatter, not by impute)
    adata.obsm["spatial"] = rng.random((n_spots, 2)).astype(np.float32)
    # Raw counts layer (OmicSage convention)
    adata.layers["counts"] = adata.X.copy()
    return adata


@pytest.fixture()
def sc_adata(rng) -> ad.AnnData:
    """100 cells × 50 genes with cell type labels."""
    n_cells, n_genes = 100, 50
    X = rng.integers(0, 100, size=(n_cells, n_genes)).astype(np.float32)
    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
    obs_names  = [f"cell_{i}" for i in range(n_cells)]
    cell_types = (["TypeA"] * 40 + ["TypeB"] * 35 + ["TypeC"] * 25)
    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({"cell_type": cell_types}, index=obs_names),
        var=pd.DataFrame(index=gene_names),
    )
    adata.layers["counts"] = adata.X.copy()
    return adata


# ---------------------------------------------------------------------------
# Helper: build a mock tangram module that does just enough
# ---------------------------------------------------------------------------

def _make_mock_tg(n_spots: int, n_genes: int, gene_names: list[str]):
    """Return a minimal mock tg namespace for Tangram calls."""
    mock_tg = types.ModuleType("tangram")

    def pp_adatas(sc, st, genes=None):
        # Tangram pre-processing: no-op for tests
        pass

    def map_cells_to_space(sc, st, **kwargs):
        # Returns an AnnData-like with obs["tg_score"] and a to_df() method
        scores = np.linspace(0.05, 0.95, n_spots)
        ad_map = MagicMock()
        ad_map.obs = pd.DataFrame({"tg_score": scores})
        # to_df() returns imputed expression
        ad_map.to_df.return_value = pd.DataFrame(
            np.random.rand(n_spots, n_genes),
            columns=gene_names,
        )
        return ad_map

    def project_genes(adata_map, adata_sc):
        # Mutates adata_map to have imputed expression; already mocked above
        pass

    mock_tg.pp_adatas = pp_adatas
    mock_tg.map_cells_to_space = map_cells_to_space
    mock_tg.project_genes = project_genes

    return mock_tg


# ---------------------------------------------------------------------------
# Import helpers
# ---------------------------------------------------------------------------

def _import_module():
    """Import spatial_impute, reloading to pick up any monkey-patching."""
    if "pipeline.modules.spatial.spatial_impute" in sys.modules:
        return importlib.reload(sys.modules["pipeline.modules.spatial.spatial_impute"])
    import pipeline.modules.spatial.spatial_impute as m
    return m


def _import_report():
    if "reports.templates.spatial.spatial_impute_report" in sys.modules:
        return importlib.reload(
            sys.modules["reports.templates.spatial.spatial_impute_report"]
        )
    import reports.templates.spatial.spatial_impute_report as r
    return r


# ---------------------------------------------------------------------------
# Module tests — skip path
# ---------------------------------------------------------------------------

class TestSkipPath:
    def test_skip_when_no_reference(self, spatial_adata):
        m = _import_module()
        adata_out, prov = m.spatial_impute(spatial_adata, adata_sc=None)
        assert prov["skipped"] is True
        assert "omicsage_spatial_impute" in adata_out.uns
        assert prov["skip_reason"] != ""

    def test_skip_does_not_mutate_input(self, spatial_adata):
        m = _import_module()
        original_uns_keys = set(spatial_adata.uns.keys())
        _, _ = m.spatial_impute(spatial_adata, adata_sc=None, inplace=False)
        assert set(spatial_adata.uns.keys()) == original_uns_keys

    def test_inplace_false_returns_copy(self, spatial_adata):
        m = _import_module()
        adata_out, _ = m.spatial_impute(spatial_adata, adata_sc=None, inplace=False)
        assert adata_out is not spatial_adata

    def test_inplace_true_modifies_original(self, spatial_adata):
        m = _import_module()
        adata_out, _ = m.spatial_impute(spatial_adata, adata_sc=None, inplace=True)
        assert adata_out is spatial_adata


# ---------------------------------------------------------------------------
# Module tests — Tangram path (mocked)
# ---------------------------------------------------------------------------

class TestTangramPath:
    def _run(self, spatial_adata, sc_adata):
        """Run spatial_impute with a mocked tangram namespace.

        Tangram is imported *inside* _run_tangram with ``import tangram as tg``,
        not at module level (because it may not be installed).  The only correct
        way to intercept that import is to inject the mock into sys.modules
        *before* the function executes, combined with patching
        ``_TANGRAM_AVAILABLE`` so the guard passes.
        """
        m = _import_module()
        gene_names = list(spatial_adata.var_names)
        mock_tg = _make_mock_tg(spatial_adata.n_obs, len(gene_names), gene_names)

        with (
            patch.object(m, "_TANGRAM_AVAILABLE", True),
            patch.dict(sys.modules, {"tangram": mock_tg}),
        ):
            return m.spatial_impute(
                spatial_adata,
                adata_sc=sc_adata,
                method="tangram",
                n_top_genes=20,
            )

    def test_output_contract_obsm(self, spatial_adata, sc_adata):
        """imputed_expression must be present in obsm as a numpy array after Tangram run."""
        adata_out, prov = self._run(spatial_adata, sc_adata)
        assert "imputed_expression" in adata_out.obsm
        assert prov["skipped"] is False
        import numpy as np
        assert isinstance(adata_out.obsm["imputed_expression"], np.ndarray)

    def test_output_shape(self, spatial_adata, sc_adata):
        """imputed_expression must be a float32 numpy array, spots as rows."""
        adata_out, _ = self._run(spatial_adata, sc_adata)
        imp = adata_out.obsm["imputed_expression"]
        import numpy as np
        assert isinstance(imp, np.ndarray)
        assert imp.shape[0] == spatial_adata.n_obs
        assert imp.dtype == np.float32

    def test_provenance_keys(self, spatial_adata, sc_adata):
        """Provenance dict must contain required keys."""
        _, prov = self._run(spatial_adata, sc_adata)
        required = {"module", "timestamp", "method", "skipped", "outputs"}
        assert required.issubset(set(prov.keys()))
        assert prov["method"] == "tangram"
        assert prov["outputs"]["n_genes_imputed"] > 0
        assert prov["outputs"]["n_spots"] == spatial_adata.n_obs

    def test_inplace_false_no_mutation(self, spatial_adata, sc_adata):
        """Input adata must not be mutated when inplace=False."""
        m = _import_module()
        gene_names = list(spatial_adata.var_names)
        mock_tg = _make_mock_tg(spatial_adata.n_obs, len(gene_names), gene_names)
        obsm_keys_before = set(spatial_adata.obsm.keys())

        with (
            patch.object(m, "_TANGRAM_AVAILABLE", True),
            patch.dict(sys.modules, {"tangram": mock_tg}),
        ):
            m.spatial_impute(spatial_adata, adata_sc=sc_adata,
                             method="tangram", n_top_genes=20, inplace=False)

        assert set(spatial_adata.obsm.keys()) == obsm_keys_before

    def test_invalid_method_raises(self, spatial_adata, sc_adata):
        m = _import_module()
        with pytest.raises(ValueError, match="Unknown imputation method"):
            m.spatial_impute(spatial_adata, adata_sc=sc_adata, method="unknown_method")


# ---------------------------------------------------------------------------
# Module tests — gimVI ImportError guard
# ---------------------------------------------------------------------------

class TestGimVIImportError:
    def test_raises_import_error_when_scvi_absent(self, spatial_adata, sc_adata):
        m = _import_module()
        with patch.object(m, "_SCVI_AVAILABLE", False):
            with pytest.raises(ImportError, match="scvi-tools"):
                m.spatial_impute(spatial_adata, adata_sc=sc_adata, method="gimvi")

    def test_raises_import_error_when_tangram_absent(self, spatial_adata, sc_adata):
        m = _import_module()
        with patch.object(m, "_TANGRAM_AVAILABLE", False):
            with pytest.raises(ImportError, match="tangram-sc"):
                m.spatial_impute(spatial_adata, adata_sc=sc_adata, method="tangram")


# ---------------------------------------------------------------------------
# Module tests — helper functions
# ---------------------------------------------------------------------------

class TestHelpers:
    def test_ensure_counts_uses_layer(self, spatial_adata):
        """_ensure_counts_in_X should prefer layers["counts"] over X."""
        m = _import_module()
        # Store different data in X vs layer to verify preference
        adata = spatial_adata.copy()
        sentinel = np.ones((adata.n_obs, adata.n_vars), dtype=np.float32) * 99
        adata.layers["counts"] = sp.csr_matrix(sentinel)
        m._ensure_counts_in_X(adata, name="test")
        result = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        np.testing.assert_array_equal(result, sentinel)

    def test_select_overlap_hvgs_returns_list(self, spatial_adata, sc_adata):
        m = _import_module()
        genes = m._select_overlap_hvgs(sc_adata, spatial_adata, n_top=10)
        assert isinstance(genes, list)
        assert len(genes) <= 10
        assert all(g in spatial_adata.var_names for g in genes)
        assert all(g in sc_adata.var_names for g in genes)

    def test_select_overlap_hvgs_no_overlap_raises(self, spatial_adata, sc_adata):
        m = _import_module()
        # Give sc_adata completely different gene names
        sc_copy = sc_adata.copy()
        sc_copy.var_names = [f"XYZ{i}" for i in range(sc_copy.n_vars)]
        with pytest.raises(ValueError, match="No overlapping genes"):
            m._select_overlap_hvgs(sc_copy, spatial_adata, n_top=10)


# ---------------------------------------------------------------------------
# Report tests
# ---------------------------------------------------------------------------

class TestReport:
    def _adata_with_impute(self, rng) -> ad.AnnData:
        """Minimal AnnData with imputation results attached."""
        n_spots, n_genes = 15, 30
        X = rng.integers(0, 50, size=(n_spots, n_genes)).astype(np.float32)
        gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
        obs_names  = [f"spot_{i}" for i in range(n_spots)]
        adata = ad.AnnData(
            X=sp.csr_matrix(X),
            obs=pd.DataFrame(index=obs_names),
            var=pd.DataFrame(index=gene_names),
        )
        adata.obsm["spatial"] = rng.random((n_spots, 2)).astype(np.float32)
        adata.layers["counts"] = adata.X.copy()
        # Attach mapping scores
        adata.obs["tangram_mapping_score"] = rng.random(n_spots).astype(np.float32)
        # Attach imputed expression as numpy array (matches module output contract).
        # Gene names are stored in uns provenance (set below).
        imputed_arr = rng.random((n_spots, n_genes)).astype(np.float32)
        adata.obsm["imputed_expression"] = imputed_arr
        # Provenance
        adata.uns["omicsage_spatial_impute"] = {
            "module":    "spatial_impute",
            "timestamp": "2026-06-11T00:00:00",
            "method":    "tangram",
            "skipped":   False,
            "skip_reason": None,
            "outputs": {
                "n_genes_imputed":    n_genes,
                "n_spots":            n_spots,
                "mean_mapping_score": 0.52,
                "n_poor_spots":       2,
                "genes_imputed":      gene_names,
                "cell_type_key":      "cell_type",
                "device":             "cpu",
            },
        }
        return adata

    def test_report_renders_without_error(self, rng, tmp_path):
        r = _import_report()
        adata = self._adata_with_impute(rng)
        out = tmp_path / "impute_report.html"
        # Patch squidpy unavailable so spatial scatter is skipped gracefully
        with patch.object(r, "_SQUIDPY_AVAILABLE", False):
            result = r.generate_spatial_impute_report(
                adata, str(out), dataset_id="test_dataset", sc_ref_label="test_ref.h5ad"
            )
        assert out.exists()
        content = out.read_text()
        assert "Spatial Imputation Report" in content
        assert "Run Summary" in content

    def test_report_skipped_state(self, rng, tmp_path):
        """Report should render a skip message when skipped=True."""
        r = _import_report()
        n_spots, n_genes = 10, 20
        X = rng.integers(0, 10, size=(n_spots, n_genes)).astype(np.float32)
        adata = ad.AnnData(X=sp.csr_matrix(X))
        adata.uns["omicsage_spatial_impute"] = {
            "module":      "spatial_impute",
            "timestamp":   "2026-06-11T00:00:00",
            "method":      "tangram",
            "skipped":     True,
            "skip_reason": "sc_reference not provided",
            "outputs":     {},
        }
        out = tmp_path / "impute_skip_report.html"
        result = r.generate_spatial_impute_report(adata, str(out))
        assert out.exists()
        content = out.read_text()
        assert "Imputation Skipped" in content or "skipped" in content.lower()

    def test_report_stat_cards_present(self, rng, tmp_path):
        r = _import_report()
        adata = self._adata_with_impute(rng)
        out = tmp_path / "impute_stat_report.html"
        with patch.object(r, "_SQUIDPY_AVAILABLE", False):
            r.generate_spatial_impute_report(
                adata, str(out), dataset_id="heart", sc_ref_label="snrna_ref.h5ad"
            )
        content = out.read_text()
        assert "stat-card" in content
        assert "Genes imputed" in content
        assert "Mean mapping score" in content

    def test_report_validation_section_present(self, rng, tmp_path):
        r = _import_report()
        adata = self._adata_with_impute(rng)
        out = tmp_path / "impute_val_report.html"
        with patch.object(r, "_SQUIDPY_AVAILABLE", False):
            r.generate_spatial_impute_report(adata, str(out))
        content = out.read_text()
        assert "Imputation Validation" in content or "Spearman" in content

    def test_report_mapping_hist_present(self, rng, tmp_path):
        r = _import_report()
        adata = self._adata_with_impute(rng)
        out = tmp_path / "impute_hist_report.html"
        with patch.object(r, "_SQUIDPY_AVAILABLE", False):
            r.generate_spatial_impute_report(adata, str(out))
        content = out.read_text()
        assert "Mapping Score" in content

    def test_report_numpy_array_obsm_handled(self, rng, tmp_path):
        """Report must handle imputed_expression stored as numpy array (checkpoint format).

        obsm["imputed_expression"] is a float32 numpy array post-checkpoint.
        Gene names live in uns["omicsage_spatial_impute"]["outputs"]["genes_imputed"].
        The report reconstructs the DataFrame from both.
        """
        r = _import_report()
        # _adata_with_impute already stores obsm["imputed_expression"] as a
        # float32 numpy array and puts gene names in uns provenance.
        # This test verifies the report can reconstruct the DataFrame from both.
        adata = self._adata_with_impute(rng)
        assert isinstance(adata.obsm["imputed_expression"], np.ndarray), (
            "Fixture must store numpy array in obsm"
        )
        out = tmp_path / "impute_numpy_report.html"
        with patch.object(r, "_SQUIDPY_AVAILABLE", False):
            r.generate_spatial_impute_report(adata, str(out))
        assert out.exists()
        content = out.read_text()
        assert "Spatial Imputation Report" in content
