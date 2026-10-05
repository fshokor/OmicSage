"""Physical neighbour graphs must never connect separate tissue sections."""
import numpy as np
import pandas as pd
import scipy.sparse as sp


def resolve_library_key(adata, library_key=None):
    key = library_key or adata.uns.get('omicsage_spatial_ingest', {}).get('library_key')
    if key:
        if key not in adata.obs or adata.obs[key].isna().any():
            raise ValueError(f'Invalid or missing tissue-section column: {key}')
        return key
    for candidate in ('patient_region_id', 'sample_name', 'library_id'):
        if candidate in adata.obs:
            if adata.obs[candidate].isna().any():
                raise ValueError(f'Missing section IDs in {candidate}')
            return candidate
    if len(adata.uns.get('spatial', {})) > 1:
        raise ValueError('Multiple spatial libraries require an explicit library_key')
    return None


def assert_section_graph(adata, library_key=None):
    key = resolve_library_key(adata, library_key)
    if key:
        labels = adata.obs[key].astype(str).to_numpy()
        for name in ('spatial_connectivities', 'spatial_distances'):
            if name not in adata.obsp:
                raise ValueError(f'{name} missing; run spatial_reduce first')
            graph = adata.obsp[name].tocoo()
            if np.any(labels[graph.row] != labels[graph.col]):
                raise ValueError('Spatial graph crosses tissue sections; rerun spatial_reduce')
    return key


def build_section_graph(adata, n_neighbors=6, coord_type=None, library_key=None):
    import anndata as ad
    import squidpy as sq
    key = resolve_library_key(adata, library_key)
    groups = list(adata.obs.groupby(key, observed=True, sort=False).indices.values()) if key else [np.arange(adata.n_obs)]
    matrices = {name: [] for name in ('spatial_connectivities', 'spatial_distances')}
    metadata = {}
    for indices in groups:
        # Coordinate-only object avoids copying a large expression matrix.
        sub = ad.AnnData(X=sp.csr_matrix((len(indices), 0)))
        sub.obsm['spatial'] = np.asarray(adata.obsm['spatial'])[indices]
        sub.uns['spatial'] = adata.uns.get('spatial', {})
        if len(indices) < 2:
            for name in matrices:
                matrices[name].append((indices, sp.csr_matrix((len(indices), len(indices)))))
            continue
        sq.gr.spatial_neighbors(sub, n_neighs=min(n_neighbors, len(indices)-1),
                               coord_type=coord_type, key_added='spatial')
        metadata = dict(sub.uns['spatial_neighbors'])
        for name in matrices:
            matrices[name].append((indices, sub.obsp[name].tocoo()))
    for name, blocks in matrices.items():
        rows, cols, values = [], [], []
        for indices, block in blocks:
            block = block.tocoo()
            rows.extend(indices[block.row]); cols.extend(indices[block.col]); values.extend(block.data)
        adata.obsp[name] = sp.csr_matrix((values, (rows, cols)), shape=(adata.n_obs, adata.n_obs))
    metadata.update({'library_key': key or '', 'n_libraries': len(groups), 'cross_library_entries': 0,
                     'connectivities_key': 'spatial_connectivities', 'distances_key': 'spatial_distances'})
    adata.uns['spatial_neighbors'] = metadata
    assert_section_graph(adata, key)
    return key
