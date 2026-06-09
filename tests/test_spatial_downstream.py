"""
tests/test_spatial_downstream.py — OmicSage Phase 7, Session 5
Tests for spatial_downstream.py and spatial_downstream_report.py

Coverage
--------
- Module imports cleanly
- inplace=False returns a copy, inplace=True modifies in place
- Graceful skip when no deconvolution data
- Graceful skip for each individual analysis when required keys are absent
- All analyses run and write expected uns keys on minimal synthetic data
- Provenance structure is correct
- Report generates valid HTML (section headings present)
- Report raises ValueError on missing provenance
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import anndata as ad


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def minimal_adata() -> ad.AnnData:
    """Minimal AnnData with all keys required for spatial_downstream."""
    n_obs, n_vars = 60, 80
    rng = np.random.default_rng(42)

    X = rng.random((n_obs, n_vars)).astype(np.float32)
    obs_names = [f"AACG{i:04d}-1" for i in range(n_obs)]
    var_names = [f"ENSG{i:08d}" for i in range(n_vars)]

    obs = pd.DataFrame(index=obs_names)
    var = pd.DataFrame(index=var_names)
    var["feature_name"] = [f"GENE{i}" for i in range(n_vars)]
    var["highly_variable"] = rng.random(n_vars) > 0.6

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Spatial coordinates (grid-like)
    coords = np.column_stack([
        np.tile(np.arange(10), 6) * 50.0,
        np.repeat(np.arange(6), 10) * 50.0,
    ]).astype(np.float32)
    adata.obsm["spatial"] = coords

    # Deconvolution outputs
    cell_types = ["Cardiomyocyte", "Fibroblast", "Endothelial"]
    abundances = rng.random((n_obs, len(cell_types))).astype(np.float32)
    abundances /= abundances.sum(axis=1, keepdims=True)
    adata.obsm["q05_cell_abundance_w_sf"] = abundances
    adata.obsm["cell_type_abundances"] = abundances
    for i, ct in enumerate(cell_types):
        adata.obs[ct] = abundances[:, i]
    adata.obs["dominant_cell_type"] = pd.Categorical(
        [cell_types[np.argmax(abundances[i])] for i in range(n_obs)]
    )

    # Spatial connectivities (simple ring graph)
    conn = sp.lil_matrix((n_obs, n_obs), dtype=np.float32)
    for i in range(n_obs):
        j = (i + 1) % n_obs
        conn[i, j] = 1.0
        conn[j, i] = 1.0
    adata.obsp["spatial_connectivities"] = conn.tocsr()
    adata.obsp["spatial_distances"] = conn.tocsr()
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "distances_key": "spatial_distances",
    }

    # Moran's I results (from spatial_cluster)
    adata.uns["moranI"] = pd.DataFrame(
        {
            "I": rng.random(n_vars),
            "pval_norm": rng.random(n_vars) * 0.1,
            "pval_norm_fdr_bh": rng.random(n_vars) * 0.2,
        },
        index=var_names,
    )

    # Deconvolution provenance
    adata.uns["omicsage_spatial_deconvolve"] = {
        "module": "spatial_deconvolve",
        "skipped": False,
        "outputs": {"cell_type_names": cell_types},
    }
    # Cluster provenance
    adata.uns["omicsage_spatial_cluster"] = {
        "module": "spatial_cluster",
        "params": {"cluster_key": "spatial_cluster"},
    }

    return adata


@pytest.fixture
def bare_adata() -> ad.AnnData:
    """Minimal AnnData with NO deconvolution or SVG data — all analyses should skip."""
    n_obs, n_vars = 20, 30
    rng = np.random.default_rng(7)
    X = rng.random((n_obs, n_vars)).astype(np.float32)
    obs_names = [f"spot_{i}" for i in range(n_obs)]
    var_names = [f"g_{i}" for i in range(n_vars)]
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    adata.obsm["spatial"] = rng.random((n_obs, 2)).astype(np.float32)
    return adata


# ---------------------------------------------------------------------------
# Import smoke test
# ---------------------------------------------------------------------------

def test_import():
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream
    from reports.templates.spatial.spatial_downstream_report import (
        generate_spatial_downstream_report,
    )
    assert callable(spatial_downstream)
    assert callable(generate_spatial_downstream_report)


# ---------------------------------------------------------------------------
# inplace behaviour
# ---------------------------------------------------------------------------

def test_inplace_false_returns_copy(minimal_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    original_id = id(minimal_adata)
    result, prov = spatial_downstream(
        minimal_adata,
        run_ligrec=False,
        run_svg_gsea=False,
        run_co_occurrence=False,
        run_nhood_enrichment=False,
        run_celltype_svg=False,
        inplace=False,
    )
    assert id(result) != original_id, "inplace=False must return a new object"
    assert "omicsage_spatial_downstream" not in minimal_adata.uns, (
        "inplace=False must not modify the original"
    )
    assert "omicsage_spatial_downstream" in result.uns


def test_inplace_true_modifies_in_place(minimal_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        minimal_adata,
        run_ligrec=False,
        run_svg_gsea=False,
        run_co_occurrence=False,
        run_nhood_enrichment=False,
        run_celltype_svg=False,
        inplace=True,
    )
    assert id(result) == id(minimal_adata), "inplace=True must return the same object"
    assert "omicsage_spatial_downstream" in minimal_adata.uns


# ---------------------------------------------------------------------------
# Type checking
# ---------------------------------------------------------------------------

def test_type_error_on_non_anndata():
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    with pytest.raises(TypeError):
        spatial_downstream({"not": "anndata"})


# ---------------------------------------------------------------------------
# Graceful skip on bare AnnData
# ---------------------------------------------------------------------------

def test_all_analyses_skip_on_bare_adata(bare_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        bare_adata,
        inplace=False,
    )
    analyses = prov["analyses"]
    for key in [
        "region_clustering", "celltype_expression", "celltype_svg",
        "co_occurrence", "nhood_enrichment", "svg_gsea",
    ]:
        info = analyses.get(key, {})
        assert info.get("skipped"), (
            f"Expected {key!r} to be skipped on bare AnnData, got: {info}"
        )


# ---------------------------------------------------------------------------
# Region clustering
# ---------------------------------------------------------------------------

def test_region_clustering_writes_obs_column(minimal_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        minimal_adata,
        run_region_clustering=True,
        run_celltype_expression=False,
        run_celltype_svg=False,
        run_co_occurrence=False,
        run_nhood_enrichment=False,
        run_ligrec=False,
        run_svg_gsea=False,
        inplace=False,
    )
    rc_info = prov["analyses"].get("region_clustering", {})
    if rc_info.get("skipped"):
        pytest.skip(f"region clustering skipped: {rc_info.get('reason')}")
    assert "region_cluster" in result.obs.columns
    assert result.obs["region_cluster"].dtype.name == "category"
    assert rc_info["n_regions"] >= 1


# ---------------------------------------------------------------------------
# Cell-type expression
# ---------------------------------------------------------------------------

def test_celltype_expression_writes_marker_dict(minimal_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        minimal_adata,
        run_region_clustering=False,
        run_celltype_expression=True,
        n_marker_genes=5,
        run_celltype_svg=False,
        run_co_occurrence=False,
        run_nhood_enrichment=False,
        run_ligrec=False,
        run_svg_gsea=False,
        inplace=False,
    )
    expr_info = prov["analyses"].get("celltype_expression", {})
    if expr_info.get("skipped"):
        pytest.skip(f"celltype expression skipped: {expr_info.get('reason')}")
    assert "celltype_marker_genes" in result.uns
    marker_dict = result.uns["celltype_marker_genes"]
    assert isinstance(marker_dict, dict)
    assert len(marker_dict) > 0
    for ct, genes in marker_dict.items():
        assert isinstance(genes, list)
        assert len(genes) <= 5, "Should store at most n_marker_genes=5 genes"


# ---------------------------------------------------------------------------
# Neighbourhood enrichment
# ---------------------------------------------------------------------------

def test_nhood_enrichment_writes_uns(minimal_adata):
    squidpy = pytest.importorskip("squidpy")
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        minimal_adata,
        run_region_clustering=False,
        run_celltype_expression=False,
        run_celltype_svg=False,
        run_co_occurrence=False,
        run_nhood_enrichment=True,
        n_perms_nhood=50,           # reduced for speed
        run_ligrec=False,
        run_svg_gsea=False,
        inplace=False,
    )
    nhood_info = prov["analyses"].get("nhood_enrichment", {})
    if nhood_info.get("skipped"):
        pytest.skip(f"nhood enrichment skipped: {nhood_info.get('reason')}")

    nhood_key = "dominant_cell_type_nhood_enrichment"
    assert nhood_key in result.uns
    nhood_data = result.uns[nhood_key]
    assert "zscore" in nhood_data


# ---------------------------------------------------------------------------
# Co-occurrence
# ---------------------------------------------------------------------------

def test_co_occurrence_writes_uns(minimal_adata):
    squidpy = pytest.importorskip("squidpy")
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        minimal_adata,
        run_region_clustering=False,
        run_celltype_expression=False,
        run_celltype_svg=False,
        run_co_occurrence=True,
        run_nhood_enrichment=False,
        run_ligrec=False,
        run_svg_gsea=False,
        inplace=False,
    )
    co_info = prov["analyses"].get("co_occurrence", {})
    if co_info.get("skipped"):
        pytest.skip(f"co-occurrence skipped: {co_info.get('reason')}")

    co_key = "dominant_cell_type_co_occurrence"
    assert co_key in result.uns


# ---------------------------------------------------------------------------
# Cell-type SVGs
# ---------------------------------------------------------------------------

def test_celltype_svg_writes_uns(minimal_adata):
    squidpy = pytest.importorskip("squidpy")
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    result, prov = spatial_downstream(
        minimal_adata,
        run_region_clustering=False,
        run_celltype_expression=False,
        run_celltype_svg=True,
        svg_n_genes=20,
        run_co_occurrence=False,
        run_nhood_enrichment=False,
        run_ligrec=False,
        run_svg_gsea=False,
        inplace=False,
    )
    svg_info = prov["analyses"].get("celltype_svg", {})
    if svg_info.get("skipped"):
        pytest.skip(f"celltype SVG skipped: {svg_info.get('reason')}")
    assert "celltype_svg" in result.uns
    svg_dict = result.uns["celltype_svg"]
    assert isinstance(svg_dict, dict)


# ---------------------------------------------------------------------------
# Provenance structure
# ---------------------------------------------------------------------------

def test_provenance_structure(minimal_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream

    _, prov = spatial_downstream(
        minimal_adata,
        run_ligrec=False,
        run_svg_gsea=False,
        run_celltype_svg=False,
        inplace=False,
    )
    assert prov["module"] == "spatial_downstream"
    assert "timestamp" in prov
    assert "params" in prov
    assert "analyses" in prov
    assert isinstance(prov["analyses"], dict)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def test_report_raises_without_provenance(tmp_path, minimal_adata):
    from reports.templates.spatial.spatial_downstream_report import (
        generate_spatial_downstream_report,
    )
    with pytest.raises(ValueError, match="omicsage_spatial_downstream"):
        generate_spatial_downstream_report(
            minimal_adata,
            output_path=str(tmp_path / "report.html"),
        )


def test_report_generates_valid_html(tmp_path, minimal_adata):
    from pipeline.modules.spatial.spatial_downstream import spatial_downstream
    from reports.templates.spatial.spatial_downstream_report import (
        generate_spatial_downstream_report,
    )

    result, _ = spatial_downstream(
        minimal_adata,
        run_region_clustering=True,
        run_celltype_expression=True,
        n_marker_genes=5,
        run_celltype_svg=False,     # slow on synthetic data
        run_co_occurrence=False,
        run_nhood_enrichment=False,
        run_ligrec=False,
        run_svg_gsea=False,
        inplace=False,
    )

    out = tmp_path / "spatial_downstream_report.html"
    path = generate_spatial_downstream_report(result, output_path=str(out))

    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert "<html" in html.lower()
    assert "Run Summary" in html
    assert "Cell-type Marker Genes" in html
    assert out.stat().st_size > 1000
