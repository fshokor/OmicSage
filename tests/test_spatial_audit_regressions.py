import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
ad=pytest.importorskip('anndata')
pytest.importorskip('squidpy')
from pipeline.modules.scripts.spatial.spatial_graph import build_section_graph, assert_section_graph
from pipeline.modules.scripts.spatial.spatial_downstream import _run_celltype_expression, _run_co_occurrence
from pipeline.modules.scripts.spatial import checkpoint_cache as cache
from pipeline.modules.scripts.spatial.spatial_qc import spatial_qc
from pipeline.modules.scripts.spatial.spatial_ingest import _load_h5ad

def two_sections():
 coords=np.array([(x,y) for x in range(3) for y in range(3)],dtype=float)
 a=ad.AnnData(sp.csr_matrix(np.ones((18,4))))
 a.obsm['spatial']=np.vstack([coords,coords])
 # Interleave sections to exercise index restoration.
 a.obs['patient_region_id']=['A']*9+['B']*9
 a.obs['dominant_cell_type']=pd.Categorical(['T','B','T']*6)
 return a[np.ravel(np.column_stack([np.arange(9),np.arange(9,18)]))].copy()

def test_overlapping_section_coordinates_never_connect():
 a=two_sections(); build_section_graph(a, n_neighbors=3, coord_type='generic')
 assert_section_graph(a)
 coo=a.obsp['spatial_connectivities'].tocoo(); labels=a.obs['patient_region_id'].to_numpy()
 assert coo.nnz>0
 assert np.all(labels[coo.row]==labels[coo.col])
 a.obsp['spatial_connectivities'][0,1]=1
 with pytest.raises(ValueError,match='crosses'):assert_section_graph(a)

def test_cooccurrence_separate_coordinate_systems(monkeypatch):
 from pipeline.modules.scripts.spatial import spatial_downstream as module
 seen=[]
 def fake(a,**kwargs):
  seen.append(a.obs['patient_region_id'].nunique())
  a.uns['dominant_cell_type_co_occurrence']={'occ':np.ones((2,2,2)),'interval':np.array([0,1,2])}
 monkeypatch.setattr(module.sq.gr,'co_occurrence',fake)
 a=two_sections();result=_run_co_occurrence(a,'dominant_cell_type',None,1)
 assert seen==[1,1];assert result['n_libraries']==2
 assert 'dominant_cell_type_co_occurrence' not in a.uns

def test_markers_exclude_negative_and_constant_correlations():
 a=ad.AnnData(np.array([[i,9-i,1] for i in range(10)],dtype=float));a.var_names=['positive','negative','constant']
 a.obs['type']=np.arange(10);a.obs['absent']=0
 a.obsm['q05_cell_abundance_w_sf']=np.column_stack([np.arange(10),np.zeros(10)])
 _run_celltype_expression(a,['type','absent'],20)
 assert a.uns['celltype_marker_genes']=={'type':['positive']}

def test_ensembl_features_use_symbols_for_mt_qc():
 a=ad.AnnData(sp.csr_matrix(np.array([[90,10],[50,50]],dtype=float)))
 a.var_names=['ENSG1','ENSG2'];a.var['feature_name']=['GENE','MT-CO1'];a.obsm['spatial']=np.zeros((2,2))
 out,meta=spatial_qc(a,min_counts=0,min_genes=0,max_mt_pct=20)
 assert out.n_obs==1;assert meta['outputs']['mt_qc_available']
 assert out.obs['pct_counts_mt'].iloc[0]==10

def test_missing_mt_report_is_explicit():
 from reports.templates.spatial.spatial_qc_report import _section_summary
 a=ad.AnnData(np.ones((1,1)))
 html=_section_summary(a,{'outputs':{'n_mt_genes_detected':0}},'test','now')
 assert 'MT QC unavailable' in html

def test_cache_rejects_changed_source_parameters_and_output(tmp_path):
 source=tmp_path/'source.h5ad';source.write_bytes(b'raw')
 cfg={'spatial':{'source':str(source),'qc':{'min_counts':500}}}
 ing=tmp_path/cache.OUTPUTS['ingest'];ing.write_bytes(b'ingested')
 qc=tmp_path/cache.OUTPUTS['qc'];qc.write_bytes(b'qc')
 cache.record('ingest',None,cfg,ing);cache.record('qc',ing,cfg,qc)
 assert cache.valid('qc',ing,cfg,qc)
 cfg['spatial']['qc']['min_counts']=600
 assert not cache.valid('qc',ing,cfg,qc)
 cfg['spatial']['qc']['min_counts']=500
 source.write_bytes(b'newraw');assert not cache.valid('qc',ing,cfg,qc)
 cache.record('ingest',None,cfg,ing);cache.record('qc',ing,cfg,qc)
 qc.write_bytes(b'corrupt');assert not cache.valid('qc',ing,cfg,qc)

def test_legacy_checkpoint_requires_rerun(tmp_path):
 output=tmp_path/cache.OUTPUTS['ingest'];output.write_bytes(b'legacy')
 assert not cache.valid('ingest',None,{},output)

def test_ligrec_qvalues_exclude_nonadjacent_section_groups(monkeypatch):
 from pipeline.modules.scripts.spatial import spatial_downstream as module
 a=two_sections(); a.obs['dominant_cell_type']=pd.Categorical(np.where(a.obs['patient_region_id']=='A','T','B'))
 build_section_graph(a,n_neighbors=3,coord_type='generic')
 def fake(data,**kwargs):
  columns=pd.MultiIndex.from_tuples([('T','T'),('T','B'),('B','T'),('B','B')])
  index=pd.MultiIndex.from_tuples([('L','R')])
  data.uns['dominant_cell_type_ligrec']={'pvalues':pd.DataFrame([[0.01]*4],index=index,columns=columns),
                                      'means':pd.DataFrame([[1.0]*4],index=index,columns=columns)}
 monkeypatch.setattr(module.sq.gr,'ligrec',fake)
 result=module._run_ligrec(a,'dominant_cell_type',10,'human',1)
 assert not result['skipped'],result
 stored=a.uns['dominant_cell_type_ligrec']
 from reports.templates.spatial.spatial_downstream_report import _deserialize_ligrec
 values=_deserialize_ligrec(stored)['qvalues'].to_numpy(dtype=float)
 assert np.isnan(values[0,1:3]).all()
 assert np.isfinite(values[0,[0,3]]).all()

def test_optional_mt_filter_still_records_percentages():
 a=ad.AnnData(sp.csr_matrix(np.array([[10,90]],dtype=float)))
 a.var_names=['GENE','MT-CO1'];a.obsm['spatial']=np.zeros((1,2))
 out,meta=spatial_qc(a,min_counts=0,min_genes=0,max_mt_pct=None)
 assert out.n_obs==1 and out.obs['pct_counts_mt'].iloc[0]==90
 assert not meta['outputs']['mt_filter_applied']
