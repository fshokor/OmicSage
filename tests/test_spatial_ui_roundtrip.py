import copy
import yaml
from ui.config_builder import build_spatial
from ui.config_io import parse_config_into_state

def test_loaded_spatial_config_keeps_all_parameters_and_nulls():
 original={'dataset_id':'heart','spatial':{
  'source':'input.h5ad','spatial_type':'h5ad','load_images':True,
  'ingest':{'library_key':'patient_region_id'},
  'qc':{'max_mt_pct':None},
  'reduce':{'coord_type':None,'normalize_total':False,'log1p':False,'flavor':'seurat'},
  'cluster':{'run_svg':False,'svg_n_genes':None,'annotation_map':{'0':'Region'}},
  'deconvolve':{'method':'nnls','library_key':None,'per_sample':True,'batch_size_st':None,'max_epochs_st':25},
  'downstream':{'n_perms_nhood':200,'ligrec_n_perms':250,'svg_n_genes':None,'region_resolution':0.7},
  'impute':{'enabled':True,'cell_type_key':'cell_type_original','max_cells_per_type':150}}}
 parsed=parse_config_into_state(original,'config.yaml')
 before=copy.deepcopy(parsed['step_params'])
 result=build_spatial('heart',parsed['data_path'],parsed['rna_path'],'Human',parsed['selected_steps'],parsed['step_params'],dataset_id='heart')
 restored=yaml.safe_load(yaml.safe_dump(result))
 for step,params in before.items():
  for key,value in params.items():
   if step=="ingest" and key in ("spatial_type","load_images"):
    assert restored["spatial"][key]==value
   else:assert restored["spatial"][step][key]==value,(step,key)
 assert before==parsed['step_params']

def test_widget_edits_override_imported_values():
 result=build_spatial('heart','input.h5ad','ref.h5ad','Human',['qc','reduce','cluster'],{
  'qc':{'max_mt_pct':None},'reduce':{'coord_type':'generic','n_comps':25},
  'cluster':{'run_svg':False,'svg_n_genes':1000,'annotation_map':{'1':'Region'}}})
 assert result['spatial']['reduce']['n_comps']==25
 assert result['spatial']['reduce']['coord_type']=='generic'
 assert result['spatial']['cluster']['run_svg'] is False
 assert result['spatial']['cluster']['annotation_map']=={'1':'Region'}


def test_real_spatial_widgets_preserve_ingestion_and_reduce_config():
 from streamlit.testing.v1 import AppTest
 import json
 original={"dataset_id":"heart","spatial":{
  "source":"input.h5ad","spatial_type":"h5ad","load_images":False,
  "ingest":{"library_key":"patient_region_id"},
  "reduce":{"n_top_genes":3000,"n_comps":50,"n_neighbors":6,
            "coord_type":None,"normalize_total":True,"target_sum":10000,
            "log1p":True,"flavor":"seurat"}}}
 parsed=parse_config_into_state(original,"config.yaml")
 at=AppTest.from_string("""
import streamlit as st
from ui._pages.p2_configure import _render_step
from ui.config_builder import build_spatial
params=st.session_state["params"]
for step in ("ingest","reduce"):
 params[step]=_render_step(step,params[step],"Spatial","Human")
st.json(build_spatial("heart","input.h5ad","","Human",["ingest","reduce"],params,dataset_id="heart"))
""")
 at.session_state["params"]=parsed["step_params"]
 at.run()
 assert not at.exception
 exported=json.loads(at.json[0].value)["spatial"]
 for key in ("spatial_type","load_images","ingest","reduce"):
  assert exported[key]==original["spatial"][key],key

def test_old_session_widget_fields_are_removed_from_export():
 cfg=build_spatial("heart","input.h5ad","","Human",["ingest","reduce"],{
  "ingest":{"library_key":"patient_region_id","spatial_type":"h5ad","load_images":True},
  "reduce":{"n_pcs_method":"elbow"}})
 assert cfg["spatial"]["ingest"]=={"library_key":"patient_region_id"}
 assert "n_pcs_method" not in cfg["spatial"]["reduce"]
