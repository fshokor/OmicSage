"""One-section selection and cell2location preparation/dispatch regressions."""
import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp
from pipeline.modules.scripts.spatial import spatial_ingest as ingest_module
from pipeline.modules.scripts.spatial import spatial_deconvolve as deconv_module


def spatial_data():
    data = ad.AnnData(sp.csr_matrix(np.ones((4, 3))),
                      obs=pd.DataFrame({'library_id': ['P1', 'P1', 'P2', 'P2']},
                                       index=['a', 'b', 'c', 'd']),
                      var=pd.DataFrame({'feature_name': ['MT-CO1', 'GENE1', 'GENE2']},
                                       index=['mt', 'g1', 'g2']))
    data.obsm['spatial'] = np.arange(8).reshape(4, 2)
    data.uns['spatial'] = {'P1': {'images': {}}, 'P2': {'images': {}}}
    data.layers['counts'] = data.X.copy()
    return data


def test_select_sample_and_image(monkeypatch):
    data = spatial_data()
    monkeypatch.setitem(ingest_module._LOADER_REGISTRY, 'h5ad',
                        (lambda *args, **kwargs: (data, 'custom', 'fixture'), True))
    out, provenance = ingest_module.spatial_ingest('fixture', spatial_type='h5ad',
                                                  library_key='library_id', sample_id='P1')
    assert out.obs_names.tolist() == ['a', 'b']
    assert list(out.uns['spatial']) == ['P1']
    assert provenance['n_obs_source'] == 4
    assert provenance['sample_id'] == 'P1'
    assert data.n_obs == 4 and len(data.uns['spatial']) == 2
    with pytest.raises(ValueError, match='Unknown sample_id'):
        ingest_module.spatial_ingest('fixture', spatial_type='h5ad',
                                     library_key='library_id', sample_id='absent')


def test_reference_overlap_excludes_mt_and_unused_layers():
    spatial = spatial_data()
    ref = ad.AnnData(sp.csr_matrix(np.arange(8).reshape(2, 4)),
                     var=pd.DataFrame(index=['g2', 'mt', 'unshared', 'g1']))
    ref.layers['counts'] = ref.X.copy()
    ref.layers['normalized_counts'] = ref.X.copy()
    prepared = deconv_module._prepare_c2l_reference(spatial, ref, 'counts')
    assert prepared.var_names.tolist() == ['g2', 'g1']
    assert list(prepared.layers) == ['counts']
    np.testing.assert_array_equal(prepared.X.toarray(), ref.layers['counts'][:, [0, 3]].toarray())


def test_joint_c2l_dispatch_and_histories(monkeypatch, tmp_path):
    spatial = spatial_data()
    ref = spatial.copy()
    history = {'elbo_train': np.array([10., 5.])}
    def fit_reference(prepared, *args, **kwargs):
        assert prepared.var_names.tolist() == ['g1', 'g2']
        prepared.uns['cell2location_reference_history'] = history
        return pd.DataFrame({'type': [1., 2.]}, index=prepared.var_names), prepared
    def fit_spatial(prepared, signatures, **kwargs):
        assert kwargs['log_every_n_epochs'] == 7
        prepared.uns['cell2location_spatial_history'] = history
        return np.ones((prepared.n_obs, 1)), ['type']
    monkeypatch.setattr(deconv_module, '_check_c2l', lambda: None)
    monkeypatch.setattr(deconv_module, '_fit_c2l_reference', fit_reference)
    monkeypatch.setattr(deconv_module, '_fit_c2l_spatial', fit_spatial)
    ref.obs['cell_type_original'] = 'type'
    out, provenance = deconv_module.spatial_deconvolve(
        spatial, ref, method='cell2location', per_sample=False, log_every_n_epochs=7)
    assert provenance['outputs']['n_shared_genes'] == 2
    assert 'cell2location_spatial_history' in out.uns
    assert 'cell2location_reference_signatures' in out.uns
    checkpoint = tmp_path / 'c2l.h5ad'
    out.write_h5ad(checkpoint)
    restored = ad.read_h5ad(checkpoint)
    np.testing.assert_array_equal(restored.uns['cell2location_spatial_history']['elbo_train'],
                                  history['elbo_train'])


def test_region_neighbors_not_limited_by_celltype_count(monkeypatch):
    from pipeline.modules.scripts.spatial import spatial_downstream as downstream
    data = ad.AnnData(sp.csr_matrix(np.ones((20, 2))))
    data.obsm['q05_cell_abundance_w_sf'] = np.ones((20, 3))
    calls = {}
    def neighbors(*args, **kwargs):
        calls.update(kwargs)
    def leiden(adata, **kwargs):
        adata.obs[kwargs['key_added']] = '0'
    def umap(adata, **kwargs):
        adata.obsm['X_umap'] = np.zeros((adata.n_obs, 2))
    monkeypatch.setattr(downstream.sc.pp, 'neighbors', neighbors)
    monkeypatch.setattr(downstream.sc.tl, 'leiden', leiden)
    monkeypatch.setattr(downstream.sc.tl, 'umap', umap)
    result = downstream._run_region_clustering(data, resolution=0.5, n_neighbors=15)
    assert not result['skipped']
    assert calls['n_neighbors'] == 15
