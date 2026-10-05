from pathlib import Path
import pytest
from pipeline.modules.scripts.spatial import checkpoint_cache as cache

def test_interrupted_write_preserves_existing_checkpoint(tmp_path):
 output=tmp_path/'checkpoint.h5ad';output.write_bytes(b'previous valid checkpoint')
 class Interrupted:
  def write_h5ad(self,path):
   Path(path).write_bytes(b'incomplete')
   raise RuntimeError('simulated interrupted serialization')
 with pytest.raises(RuntimeError,match='interrupted'):cache.atomic_write_h5ad(Interrupted(),output)
 assert output.read_bytes()==b'previous valid checkpoint'
 assert list(tmp_path.iterdir())==[output]

def test_completed_write_replaces_checkpoint(tmp_path):
 output=tmp_path/'checkpoint.h5ad';output.write_bytes(b'old')
 class Complete:
  def write_h5ad(self,path):Path(path).write_bytes(b'new complete checkpoint')
 cache.atomic_write_h5ad(Complete(),output)
 assert output.read_bytes()==b'new complete checkpoint'
 assert list(tmp_path.iterdir())==[output]

def test_error_identifies_corrupt_ancestor(tmp_path):
 cfg={};ingest=tmp_path/cache.OUTPUTS['ingest'];ingest.write_bytes(b'good')
 cache.record('ingest',None,cfg,ingest)
 ingest.write_bytes(b'corrupt')
 target=tmp_path/cache.OUTPUTS['downstream']
 with pytest.raises(ValueError,match="checksum changed.*01_ingested.*Rebuild from 'ingest'"):
  cache.require_predecessor('downstream',cfg,target)
