import json
import yaml
from streamlit.testing.v1 import AppTest
from ui.config_builder import build_spatial

def test_disabled_mt_survives_yaml_export():
 cfg=build_spatial('heart','input.h5ad','ref.h5ad','Human',['ingest','qc'],{'qc':{'max_mt_pct':None}})
 restored=yaml.safe_load(yaml.safe_dump(cfg))
 assert 'max_mt_pct' in restored['spatial']['qc']
 assert restored['spatial']['qc']['max_mt_pct'] is None

def test_numeric_mt_survives_yaml_export():
 cfg=build_spatial('heart','input.h5ad','ref.h5ad','Human',['ingest','qc'],{'qc':{'max_mt_pct':35.5}})
 assert cfg['spatial']['qc']['max_mt_pct']==35.5

def test_spatial_qc_page_with_disabled_mt_and_toggle():
 at=AppTest.from_string('''
import streamlit as st
from ui._pages.p2_configure import _render_step
params=_render_step("qc", {"max_mt_pct": None}, "Spatial", "Human")
st.json(params)
''').run()
 assert not at.exception
 assert json.loads(at.json[0].value)['max_mt_pct'] is None
 at.toggle[0].set_value(True).run()
 assert not at.exception
 assert json.loads(at.json[0].value)['max_mt_pct']==20.0
 at.toggle[0].set_value(False).run()
 assert not at.exception
 assert json.loads(at.json[0].value)['max_mt_pct'] is None
