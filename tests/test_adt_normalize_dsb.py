"""
tests/test_adt_normalize_dsb.py

Tests for the DSB normalization path added to normalize_adt().

Run with:
    python -m pytest tests/test_adt_normalize_dsb.py -v --tb=short
"""

from __future__ import annotations

import numpy as np
import pytest

try:
    import muon  # noqa: F401
    MUON_AVAILABLE = True
except ImportError:
    MUON_AVAILABLE = False

from anndata import AnnData


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_adata(
    n_cells: int = 80,
    n_proteins: int = 20,
    seed: int = 0,
    protein_names: list[str] | None = None,
) -> AnnData:
    """
    Synthetic ADT AnnData with raw integer counts.
    Uses two well-separated Gaussian clusters so Leiden produces ≥ 2 groups.
    """
    rng = np.random.default_rng(seed)
    half = n_cells // 2
    # Group A: higher counts for first 10 proteins
    group_a = rng.poisson(lam=30, size=(half, n_proteins)).astype(np.float32)
    group_a[:, :10] += rng.poisson(lam=50, size=(half, 10))
    # Group B: higher counts for last 10 proteins
    group_b = rng.poisson(lam=30, size=(n_cells - half, n_proteins)).astype(np.float32)
    group_b[:, 10:] += rng.poisson(lam=50, size=(n_cells - half, 10))

    X = np.vstack([group_a, group_b])
    if protein_names is None:
        protein_names = [f"CD{i}" for i in range(n_proteins)]
    adata = AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = protein_names
    return adata


def _make_empty_adata(
    n_empty: int = 200,
    n_proteins: int = 20,
    seed: int = 42,
    protein_names: list[str] | None = None,
) -> AnnData:
    """
    Synthetic empty-droplet AnnData with low Poisson counts (ambient signal only).
    """
    rng = np.random.default_rng(seed)
    X = rng.poisson(lam=3, size=(n_empty, n_proteins)).astype(np.float32)
    if protein_names is None:
        protein_names = [f"CD{i}" for i in range(n_proteins)]
    empty = AnnData(X=X)
    empty.obs_names = [f"empty_{i}" for i in range(n_empty)]
    empty.var_names = protein_names
    return empty


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def adata_raw():
    return _make_adata()


@pytest.fixture
def empty_adata():
    return _make_empty_adata()


# ---------------------------------------------------------------------------
# Validation tests (existing contract — no regressions)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not MUON_AVAILABLE, reason="muon not installed")
class TestNormalizeAdtValidation:

    def test_rejects_non_anndata(self):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        with pytest.raises(TypeError, match="AnnData"):
            normalize_adt("not_an_anndata")

    def test_rejects_empty_cells(self):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        empty = AnnData(X=np.zeros((0, 5), dtype=float))
        empty.var_names = [f"CD{i}" for i in range(5)]
        with pytest.raises(ValueError, match="0 cells"):
            normalize_adt(empty)

    def test_rejects_empty_proteins(self):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        no_prot = AnnData(X=np.zeros((10, 0), dtype=float))
        with pytest.raises(ValueError, match="0 proteins"):
            normalize_adt(no_prot)

    def test_rejects_bad_clr_axis(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        with pytest.raises(ValueError, match="clr_axis"):
            normalize_adt(adata_raw, clr_axis=2)

    def test_rejects_dsb_empty_non_anndata(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        with pytest.raises(TypeError, match="AnnData"):
            normalize_adt(adata_raw, dsb_empty_adata="bad")

    def test_rejects_dsb_empty_zero_rows(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        zero_empty = AnnData(X=np.zeros((0, adata_raw.n_vars), dtype=float))
        zero_empty.var_names = list(adata_raw.var_names)
        with pytest.raises(ValueError, match="0 rows"):
            normalize_adt(adata_raw, dsb_empty_adata=zero_empty)

    def test_rejects_dsb_protein_mismatch(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        bad_empty = _make_empty_adata(n_proteins=adata_raw.n_vars + 1)
        with pytest.raises(ValueError, match="proteins"):
            normalize_adt(adata_raw, dsb_empty_adata=bad_empty)

    def test_rejects_dsb_varnames_mismatch(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        bad_empty = _make_empty_adata(
            n_proteins=adata_raw.n_vars,
            protein_names=[f"OTHER{i}" for i in range(adata_raw.n_vars)],
        )
        with pytest.raises(ValueError, match="var_names"):
            normalize_adt(adata_raw, dsb_empty_adata=bad_empty)


# ---------------------------------------------------------------------------
# CLR-only tests (regression — no DSB)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not MUON_AVAILABLE, reason="muon not installed")
class TestNormalizeAdtClrOnly:

    def test_clr_layers_written(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, metrics = normalize_adt(adata_raw)
        assert "counts" in out.layers
        assert "adt_clr" in out.layers
        assert "adt_dsb" not in out.layers

    def test_clr_x_equals_adt_clr(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, _ = normalize_adt(adata_raw)
        np.testing.assert_array_almost_equal(
            np.asarray(out.X), np.asarray(out.layers["adt_clr"])
        )

    def test_clr_counts_preserved(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        raw_X = adata_raw.X.copy()
        out, _ = normalize_adt(adata_raw)
        np.testing.assert_array_equal(np.asarray(out.layers["counts"]), raw_X)

    def test_clr_metrics_dsb_false(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        _, metrics = normalize_adt(adata_raw)
        assert metrics["dsb_applied"] is False
        assert metrics["active_layer"] == "adt_clr"

    def test_inplace_false_does_not_mutate(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        raw_X = adata_raw.X.copy()
        normalize_adt(adata_raw, inplace=False)
        np.testing.assert_array_equal(np.asarray(adata_raw.X), raw_X)

    def test_inplace_true_mutates(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        raw_X = adata_raw.X.copy()
        normalize_adt(adata_raw, inplace=True)
        # After CLR, .X should differ from raw integer counts
        assert not np.allclose(np.asarray(adata_raw.X), raw_X)

    def test_clr_axis_1(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, metrics = normalize_adt(adata_raw, clr_axis=1)
        assert metrics["clr_axis"] == 1
        assert "adt_clr" in out.layers

    def test_provenance_written(self, adata_raw):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, _ = normalize_adt(adata_raw)
        uns = out.uns["omicsage_adt_normalize"]
        assert "timestamp" in uns
        assert uns["dsb_applied"] is False
        assert uns["active_layer"] == "adt_clr"


# ---------------------------------------------------------------------------
# DSB tests (new)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not MUON_AVAILABLE, reason="muon not installed")
class TestNormalizeAdtDsb:

    def test_dsb_layers_written(self, adata_raw, empty_adata):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, metrics = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        assert "adt_dsb" in out.layers, "adt_dsb layer must be written when DSB runs"
        assert "adt_clr" in out.layers, "adt_clr layer must still be present after DSB"
        assert "counts" in out.layers, "raw counts layer must be preserved"

    def test_dsb_x_equals_adt_dsb(self, adata_raw, empty_adata):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, _ = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        np.testing.assert_array_almost_equal(
            np.asarray(out.X), np.asarray(out.layers["adt_dsb"])
        )

    def test_dsb_raw_counts_preserved(self, adata_raw, empty_adata):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        raw_X = adata_raw.X.copy()
        out, _ = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        np.testing.assert_array_equal(np.asarray(out.layers["counts"]), raw_X)

    def test_dsb_metrics_populated(self, adata_raw, empty_adata):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        _, metrics = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        assert metrics["dsb_applied"] is True
        assert metrics["active_layer"] == "adt_dsb"
        assert metrics["dsb_n_empty_droplets"] == empty_adata.n_obs
        assert "dsb_mean" in metrics
        assert "dsb_fraction_positive" in metrics
        assert 0.0 <= metrics["dsb_fraction_positive"] <= 1.0

    def test_dsb_x_differs_from_clr(self, adata_raw, empty_adata):
        """After DSB, adata.X (DSB) must differ from adt_clr."""
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, _ = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        x_dsb = np.asarray(out.X)
        x_clr = np.asarray(out.layers["adt_clr"])
        assert not np.allclose(x_dsb, x_clr), \
            "DSB and CLR values should differ — DSB subtracts ambient background"

    def test_dsb_provenance_in_uns(self, adata_raw, empty_adata):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, _ = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        uns = out.uns["omicsage_adt_normalize"]
        assert uns["dsb_applied"] is True
        assert uns["active_layer"] == "adt_dsb"
        assert uns["normalized_in_layer"] == "adt_dsb"
        assert "dsb_params" in uns

    def test_dsb_isotype_controls_present(self, adata_raw, empty_adata):
        """When valid isotype controls are provided, they're recorded in metrics."""
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        # Use the first 2 proteins as fake isotype controls
        iso = list(adata_raw.var_names[:2])
        _, metrics = normalize_adt(
            adata_raw,
            dsb_empty_adata=empty_adata,
            isotype_controls=iso,
            dsb_denoise=False,
        )
        assert metrics["dsb_isotype_controls_used"] == iso

    def test_dsb_isotype_controls_missing_silently_filtered(
        self, adata_raw, empty_adata
    ):
        """Controls not in var_names should be silently filtered (no exception)."""
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        iso = ["NONEXISTENT_CTRL_1", "NONEXISTENT_CTRL_2"]
        _, metrics = normalize_adt(
            adata_raw,
            dsb_empty_adata=empty_adata,
            isotype_controls=iso,
            dsb_denoise=False,
        )
        # All filtered out → empty list used
        assert metrics["dsb_isotype_controls_used"] == []

    def test_dsb_no_muon_layer_name_collision(self, adata_raw, empty_adata):
        """The muon 'dsb' layer should be renamed to 'adt_dsb'."""
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        out, _ = normalize_adt(
            adata_raw, dsb_empty_adata=empty_adata, dsb_denoise=False
        )
        assert "dsb" not in out.layers, \
            "muon's raw 'dsb' layer should be renamed to 'adt_dsb'"
        assert "adt_dsb" in out.layers

    def test_dsb_inplace_false_does_not_mutate(self, adata_raw, empty_adata):
        from pipeline.modules.scripts.cite.adt_normalize import normalize_adt
        raw_X = adata_raw.X.copy()
        normalize_adt(adata_raw, dsb_empty_adata=empty_adata,
                      dsb_denoise=False, inplace=False)
        np.testing.assert_array_equal(np.asarray(adata_raw.X), raw_X)
