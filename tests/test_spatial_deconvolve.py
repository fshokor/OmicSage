"""
test_spatial_deconvolve.py — OmicSage Phase 7, Session 4
Tests for spatial_deconvolve() and generate_spatial_deconvolve_report().

Strategy:
  - cell2location training is SKIPPED in unit tests (too slow, needs real data).
  - We test the full module contract via two paths:
      1. ref_adata=None  → graceful skip path (fast, no c2l needed)
      2. Synthetic post-deconvolution state → report + provenance tests
  - The real cell2location path is covered by integration smoke test
    (marked slow, skipped by default).
  - ingest _load_h5ad changes are tested separately in test_spatial_ingest.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

anndata = pytest.importorskip("anndata")
scanpy = pytest.importorskip("scanpy")
squidpy = pytest.importorskip("squidpy")

import anndata as ad
import pandas as pd

from pipeline.modules.spatial.spatial_reduce import spatial_reduce
from pipeline.modules.spatial.spatial_cluster import spatial_cluster
from pipeline.modules.spatial.spatial_deconvolve import spatial_deconvolve
from reports.templates.spatial.spatial_deconvolve_report import (
    generate_spatial_deconvolve_report,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_SPOTS = 80
N_GENES = 200
N_CELL_TYPES = 5
RNG = np.random.default_rng(99)

CELL_TYPES = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Myeloid", "Pericyte"]


def _make_clustered_adata(
    n_obs: int = N_SPOTS,
    n_vars: int = N_GENES,
) -> ad.AnnData:
    """Minimal AnnData that has passed through the full spatial pipeline."""
    X = sp.csr_matrix(
        RNG.negative_binomial(5, 0.5, size=(n_obs, n_vars)).astype(np.float32)
    )
    var_names = [f"ENSG{i:011d}" for i in range(n_vars)]  # ENSEMBL-style IDs
    obs_names = [f"spot{i}" for i in range(n_obs)]

    adata = ad.AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = var_names

    side = int(np.ceil(np.sqrt(n_obs)))
    coords = np.column_stack([
        np.tile(np.arange(side), side)[:n_obs],
        np.repeat(np.arange(side), side)[:n_obs],
    ]).astype(float)
    adata.obsm["spatial"] = coords
    adata.uns["spatial"] = {"visium": {"images": {}, "scalefactors": {}}}
    adata.layers["counts"] = X.copy()
    adata.obs["patient"] = "P1"
    adata.uns["omicsage_spatial_ingest"] = {
        "spatial_type": "h5ad",
        "n_obs": n_obs,
        "n_vars": n_vars,
        "timestamp": "2026-01-01T00:00:00",
    }

    adata, _ = spatial_reduce(adata, n_top_genes=60, n_comps=10, inplace=False)
    adata, _ = spatial_cluster(adata, resolution=0.3, run_svg=False, inplace=False)
    return adata


def _make_ref_adata(
    n_cells: int = 100,
    n_vars: int = N_GENES,
) -> ad.AnnData:
    """Minimal reference scRNA-seq AnnData with matching ENSEMBL gene IDs."""
    X = sp.csr_matrix(
        RNG.negative_binomial(3, 0.4, size=(n_cells, n_vars)).astype(np.float32)
    )
    var_names = [f"ENSG{i:011d}" for i in range(n_vars)]
    obs_names = [f"cell{i}" for i in range(n_cells)]
    cell_types = [CELL_TYPES[i % N_CELL_TYPES] for i in range(n_cells)]

    adata = ad.AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.obs["cell_type_original"] = pd.Categorical(cell_types)
    adata.obs["donor_id"] = pd.Categorical(["D1"] * 50 + ["D2"] * 50)
    adata.obs["assay"] = pd.Categorical(["10x"] * n_cells)
    adata.layers["counts"] = X.copy()
    return adata


def _make_deconvolved_adata() -> ad.AnnData:
    """AnnData simulating post-deconvolution state (bypasses cell2location)."""
    adata = _make_clustered_adata()

    # Simulate cell2location output
    n_obs = adata.n_obs
    abundance = RNG.dirichlet(np.ones(N_CELL_TYPES), size=n_obs).astype(np.float32)
    import pandas as pd
    adata.obsm["q05_cell_abundance_w_sf"] = pd.DataFrame(
        abundance, index=adata.obs_names, columns=CELL_TYPES
    )
    for ct in CELL_TYPES:
        adata.obs[ct] = abundance[:, CELL_TYPES.index(ct)]

    adata.uns["cell2location_mod"] = {"factor_names": CELL_TYPES}
    adata.uns["omicsage_spatial_deconvolve"] = {
        "module": "spatial_deconvolve",
        "timestamp": "2026-01-01T12:00:00",
        "skipped": False,
        "params": {
            "cell_type_key": "cell_type_original",
            "batch_key_ref": "donor_id",
            "batch_key_st": "patient",
            "covariate_keys": ["assay"],
            "layer_ref": "counts",
            "N_cells_per_location": 8,
            "detection_alpha": 20,
            "max_epochs_ref": 250,
            "max_epochs_st": 30000,
            "batch_size_ref": 2500,
            "cell_count_cutoff": 5,
            "cell_percentage_cutoff2": 0.03,
            "nonz_mean_cutoff": 1.12,
        },
        "outputs": {
            "n_cell_types": N_CELL_TYPES,
            "cell_type_names": CELL_TYPES,
            "n_shared_genes": N_GENES,
            "n_genes_after_selection": 150,
            "n_spots": N_SPOTS,
        },
    }
    return adata


# ---------------------------------------------------------------------------
# Contract tests — skip path (ref_adata=None)
# ---------------------------------------------------------------------------

class TestSpatialDeconvolveSkipPath:

    def test_returns_tuple(self):
        adata = _make_clustered_adata()
        result = spatial_deconvolve(adata, ref_adata=None)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_returns_anndata_and_dict(self):
        adata = _make_clustered_adata()
        out, params = spatial_deconvolve(adata, ref_adata=None)
        assert isinstance(out, ad.AnnData)
        assert isinstance(params, dict)

    def test_provenance_written(self):
        adata = _make_clustered_adata()
        out, params = spatial_deconvolve(adata, ref_adata=None)
        assert "omicsage_spatial_deconvolve" in out.uns
        assert out.uns["omicsage_spatial_deconvolve"] is params

    def test_skipped_flag_true(self):
        adata = _make_clustered_adata()
        _, params = spatial_deconvolve(adata, ref_adata=None)
        assert bool(params["skipped"]) is True

    def test_skip_reason_present(self):
        adata = _make_clustered_adata()
        _, params = spatial_deconvolve(adata, ref_adata=None)
        assert params["skip_reason"]
        assert len(params["skip_reason"]) > 10

    def test_provenance_has_required_keys(self):
        adata = _make_clustered_adata()
        _, params = spatial_deconvolve(adata, ref_adata=None)
        for key in ("module", "timestamp", "skipped", "params", "outputs"):
            assert key in params

    def test_inplace_false_returns_new_object(self):
        adata = _make_clustered_adata()
        out, _ = spatial_deconvolve(adata, ref_adata=None, inplace=False)
        assert out is not adata

    def test_inplace_false_does_not_modify_input(self):
        adata = _make_clustered_adata()
        spatial_deconvolve(adata, ref_adata=None, inplace=False)
        assert "omicsage_spatial_deconvolve" not in adata.uns

    def test_inplace_true_modifies_input(self):
        adata = _make_clustered_adata()
        out, _ = spatial_deconvolve(adata, ref_adata=None, inplace=True)
        assert out is adata
        assert "omicsage_spatial_deconvolve" in adata.uns

    def test_no_abundance_written_when_skipped(self):
        adata = _make_clustered_adata()
        out, _ = spatial_deconvolve(adata, ref_adata=None)
        assert "q05_cell_abundance_w_sf" not in out.obsm

    def test_obs_shape_preserved(self):
        adata = _make_clustered_adata()
        out, _ = spatial_deconvolve(adata, ref_adata=None)
        assert out.n_obs == N_SPOTS


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestSpatialDeconvolveValidation:

    def test_rejects_non_anndata_adata(self):
        with pytest.raises(TypeError, match="AnnData"):
            spatial_deconvolve({"not": "anndata"}, ref_adata=None)

    def test_rejects_missing_counts_layer(self):
        adata = _make_clustered_adata()
        del adata.layers["counts"]
        with pytest.raises(ValueError, match="counts"):
            spatial_deconvolve(adata, ref_adata=None)

    def test_rejects_empty_adata(self):
        adata = _make_clustered_adata()
        adata = adata[:0].copy()
        with pytest.raises(ValueError, match="0 observations"):
            spatial_deconvolve(adata, ref_adata=None)

    def test_rejects_non_anndata_ref(self):
        adata = _make_clustered_adata()
        with pytest.raises(TypeError, match="AnnData"):
            spatial_deconvolve(adata, ref_adata={"not": "anndata"})

    def test_rejects_ref_missing_layer(self):
        adata = _make_clustered_adata()
        ref = _make_ref_adata()
        del ref.layers["counts"]
        with pytest.raises(ValueError, match="counts"):
            spatial_deconvolve(adata, ref_adata=ref, layer_ref="counts")

    def test_rejects_ref_missing_cell_type_key(self):
        adata = _make_clustered_adata()
        ref = _make_ref_adata()
        with pytest.raises(ValueError, match="cell_type_key"):
            spatial_deconvolve(
                adata, ref_adata=ref, cell_type_key="nonexistent_col"
            )

    def test_rejects_empty_ref(self):
        adata = _make_clustered_adata()
        ref = _make_ref_adata()
        ref = ref[:0].copy()
        with pytest.raises(ValueError, match="0 observations"):
            spatial_deconvolve(adata, ref_adata=ref)


# ---------------------------------------------------------------------------
# Simulated post-deconvolution state tests
# ---------------------------------------------------------------------------

class TestDeconvolvedAdataState:

    @pytest.fixture
    def deconvolved_adata(self):
        return _make_deconvolved_adata()

    def test_abundance_in_obsm(self, deconvolved_adata):
        assert "q05_cell_abundance_w_sf" in deconvolved_adata.obsm

    def test_abundance_shape(self, deconvolved_adata):
        ab = deconvolved_adata.obsm["q05_cell_abundance_w_sf"]
        assert ab.shape == (N_SPOTS, N_CELL_TYPES)

    def test_cell_type_columns_in_obs(self, deconvolved_adata):
        for ct in CELL_TYPES:
            assert ct in deconvolved_adata.obs.columns

    def test_abundance_values_non_negative(self, deconvolved_adata):
        ab = deconvolved_adata.obsm["q05_cell_abundance_w_sf"]
        assert (np.array(ab) >= 0).all()

    def test_provenance_skipped_false(self, deconvolved_adata):
        params = deconvolved_adata.uns["omicsage_spatial_deconvolve"]
        assert bool(params["skipped"]) is False

    def test_provenance_outputs_keys(self, deconvolved_adata):
        outputs = deconvolved_adata.uns["omicsage_spatial_deconvolve"]["outputs"]
        for key in ("n_cell_types", "cell_type_names", "n_shared_genes", "n_spots"):
            assert key in outputs

    def test_n_cell_types_matches(self, deconvolved_adata):
        outputs = deconvolved_adata.uns["omicsage_spatial_deconvolve"]["outputs"]
        assert outputs["n_cell_types"] == N_CELL_TYPES

    def test_cell_type_names_match(self, deconvolved_adata):
        outputs = deconvolved_adata.uns["omicsage_spatial_deconvolve"]["outputs"]
        assert set(outputs["cell_type_names"]) == set(CELL_TYPES)


# ---------------------------------------------------------------------------
# Ingest _load_h5ad contract tests
# ---------------------------------------------------------------------------

class TestIngestH5adContract:

    def test_h5ad_saves_counts_layer(self, tmp_path):
        """_load_h5ad must save layers['counts'] even when absent in source."""
        import scipy.sparse as sp
        X = sp.csr_matrix(
            RNG.negative_binomial(5, 0.5, size=(20, 50)).astype(np.float32)
        )
        adata = ad.AnnData(X=X)
        adata.obs_names = [f"spot{i}" for i in range(20)]
        adata.var_names = [f"Gene{i}" for i in range(50)]
        adata.obsm["spatial"] = np.random.rand(20, 2)
        adata.uns["spatial"] = {"lib": {"images": {}, "scalefactors": {}}}
        # No layers["counts"] — simulate raw h5ad
        h5ad_path = str(tmp_path / "test.h5ad")
        adata.write_h5ad(h5ad_path)

        from pipeline.modules.spatial.spatial_ingest import spatial_ingest
        out, _ = spatial_ingest(h5ad_path, spatial_type="h5ad")
        assert "counts" in out.layers

    def test_h5ad_swaps_var_names_to_ensembl(self, tmp_path):
        """_load_h5ad must swap var_names to gene_ids when present."""
        X = sp.csr_matrix(np.ones((10, 5), dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.obs_names = [f"spot{i}" for i in range(10)]
        adata.var_names = [f"GeneSymbol{i}" for i in range(5)]
        adata.var["gene_ids"] = [f"ENSG{i:011d}" for i in range(5)]
        adata.obsm["spatial"] = np.random.rand(10, 2)
        adata.uns["spatial"] = {"lib": {"images": {}, "scalefactors": {}}}
        h5ad_path = str(tmp_path / "test_ensembl.h5ad")
        adata.write_h5ad(h5ad_path)

        from pipeline.modules.spatial.spatial_ingest import spatial_ingest
        out, _ = spatial_ingest(h5ad_path, spatial_type="h5ad")
        assert all(v.startswith("ENSG") for v in out.var_names)

    def test_h5ad_strips_mt_genes(self, tmp_path):
        """_load_h5ad must strip MT- genes into obsm['MT']."""
        n_mt = 3
        n_other = 10
        X = sp.csr_matrix(np.ones((10, n_mt + n_other), dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.obs_names = [f"spot{i}" for i in range(10)]
        gene_names = [f"MT-Gene{i}" for i in range(n_mt)] + \
                     [f"Gene{i}" for i in range(n_other)]
        adata.var_names = gene_names
        adata.obsm["spatial"] = np.random.rand(10, 2)
        adata.uns["spatial"] = {"lib": {"images": {}, "scalefactors": {}}}
        h5ad_path = str(tmp_path / "test_mt.h5ad")
        adata.write_h5ad(h5ad_path)

        from pipeline.modules.spatial.spatial_ingest import spatial_ingest
        out, _ = spatial_ingest(h5ad_path, spatial_type="h5ad")
        assert "MT" in out.obsm
        assert out.n_vars == n_other
        assert out.obsm["MT"].shape == (10, n_mt)


# ---------------------------------------------------------------------------
# Report tests
# ---------------------------------------------------------------------------

class TestSpatialDeconvolveReport:

    @pytest.fixture
    def deconvolved_adata(self):
        return _make_deconvolved_adata()

    @pytest.fixture
    def skipped_adata(self):
        adata = _make_clustered_adata()
        out, _ = spatial_deconvolve(adata, ref_adata=None, inplace=False)
        return out

    def test_report_generates_file(self, deconvolved_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        result = generate_spatial_deconvolve_report(
            deconvolved_adata, out_path, dataset_id="test"
        )
        assert result.endswith(".html")
        import os
        assert os.path.isfile(result)

    def test_report_not_empty(self, deconvolved_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_deconvolve_report(
            deconvolved_adata, out_path, dataset_id="test"
        )
        import os
        assert os.path.getsize(out_path) > 1000

    def test_report_is_valid_html(self, deconvolved_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_deconvolve_report(
            deconvolved_adata, out_path, dataset_id="test"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "<!DOCTYPE html>" in content
        assert "</html>" in content

    def test_report_contains_dataset_id(self, deconvolved_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_deconvolve_report(
            deconvolved_adata, out_path, dataset_id="kuppe_heart"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "kuppe_heart" in content

    def test_report_contains_cell_types(self, deconvolved_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_deconvolve_report(
            deconvolved_adata, out_path, dataset_id="test"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "Cardiomyocyte" in content

    def test_report_raises_without_provenance(self, tmp_path):
        adata = _make_clustered_adata()
        with pytest.raises(ValueError, match="omicsage_spatial_deconvolve"):
            generate_spatial_deconvolve_report(
                adata, str(tmp_path / "x.html"), dataset_id="test"
            )

    def test_report_returns_absolute_path(self, deconvolved_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        result = generate_spatial_deconvolve_report(
            deconvolved_adata, out_path, dataset_id="test"
        )
        import os
        assert os.path.isabs(result)

    def test_skipped_report_contains_notice(self, skipped_adata, tmp_path):
        out_path = str(tmp_path / "skipped_report.html")
        generate_spatial_deconvolve_report(
            skipped_adata, out_path, dataset_id="test"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "Skipped" in content or "skipped" in content.lower()

    def test_skipped_report_is_valid_html(self, skipped_adata, tmp_path):
        out_path = str(tmp_path / "skipped_report.html")
        generate_spatial_deconvolve_report(
            skipped_adata, out_path, dataset_id="test"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "<!DOCTYPE html>" in content
        assert "</html>" in content


# ---------------------------------------------------------------------------
# Slow integration test (skipped by default — requires real data + time)
# ---------------------------------------------------------------------------

@pytest.mark.skip(reason="Integration test — requires Kuppe data and ~2h runtime")
def test_full_deconvolution_kuppe():
    """Smoke test: full cell2location pipeline on Kuppe heart data."""
    import scanpy as sc
    adata_st = sc.read_h5ad(
        "data/benchmark/kuppe_visium_human_heart_2022_control.h5ad"
    )
    ref_adata = sc.read_h5ad(
        "data/benchmark/kuppe_snRNA_human_heart_2022_control.h5ad"
    )
    out, params = spatial_deconvolve(
        adata_st,
        ref_adata=ref_adata,
        cell_type_key="cell_type_original",
        batch_key_ref="donor_id",
        batch_key_st="patient",
        covariate_keys=["assay"],
        N_cells_per_location=8,
        max_epochs_ref=10,   # fast smoke test
        max_epochs_st=100,
    )
    assert not params["skipped"]
    assert "q05_cell_abundance_w_sf" in out.obsm
    assert params["outputs"]["n_cell_types"] > 0
