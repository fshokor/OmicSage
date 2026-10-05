"""Content-addressed checkpoints, including recursively validated predecessors."""
import hashlib
import json
from pathlib import Path

OUTPUTS = dict(zip(('ingest','qc','reduce','cluster','deconvolve','downstream','impute'),
                  ('01_ingested.h5ad','02_qc.h5ad','03_reduced.h5ad','04_clustered.h5ad',
                   '05_deconvolved.h5ad','06_downstream.h5ad','07_imputed.h5ad')))
PARENTS = {'qc':'ingest','reduce':'qc','cluster':'reduce','deconvolve':'cluster',
           'downstream':'deconvolve','impute':'cluster'}
_hashes = {}

def digest(path):
    path = Path(path).resolve(); stat = path.stat()
    key = (str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    if key not in _hashes:
        h = hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024), b''):h.update(block)
        _hashes[key] = h.hexdigest()
    return _hashes[key]

def manifest_path(output):return Path(str(output)+'.cache.json')

def signature(step, input_path, cfg):
    spatial = cfg.get('spatial', {})
    own = ({k:v for k,v in spatial.items() if k in ('source','spatial_type','counts_file','library_id','library_key','load_images','ingest')}
           if step == 'ingest' else spatial.get(step, {}))
    paths = [input_path] if input_path else []
    if step == 'ingest':
        source = spatial.get('source')
        if source and Path(source).exists():
            source = Path(source)
            paths.extend([source] if source.is_file() else sorted(p for p in source.rglob('*') if p.is_file()))
    if isinstance(own, dict):
        for name in ('ref_path','sc_reference_path'):
            if own.get(name):paths.append(Path(own[name]))
    code = Path(__file__).resolve().parent
    root = code.parents[3]
    sources = list(code.glob('*.py')) + list((root/'reports/templates/spatial').glob('*.py'))
    runner = root/'run_spatial_pipeline.py'
    if runner.exists():sources.append(runner)
    record = {'schema':1,'step':step,'config':own,'dataset_id':cfg.get('dataset_id'),
              'files':{str(Path(p).resolve()):digest(p) for p in paths if p},
              'code':{str(p.resolve()):digest(p) for p in sorted(sources)}}
    return hashlib.sha256(json.dumps(record,sort_keys=True,default=str).encode()).hexdigest()

def predecessor(step, cfg, output):
    parent = PARENTS.get(step)
    if parent is None:return None, None
    directory = Path(output).parent
    if step == 'downstream' and (not (directory/OUTPUTS[parent]).exists()
                                or not cfg.get('spatial', {}).get('deconvolve', {}).get('enabled', True)):
        parent = 'cluster'
    return parent, directory/OUTPUTS[parent]

def valid(step, input_path, cfg, output):
    output = Path(output)
    if not output.exists() or not manifest_path(output).exists():return False
    try:
        parent, source = predecessor(step,cfg,output)
        if parent and not valid(parent,predecessor(parent,cfg,source)[1],cfg,source):return False
        saved = json.loads(manifest_path(output).read_text())
        return saved['signature']==signature(step,input_path,cfg) and saved['output_sha256']==digest(output)
    except (OSError,ValueError,KeyError,TypeError):return False

def first_invalid(step, input_path, cfg, output):
    """Report the earliest failed dependency, not just the requested parent."""
    output = Path(output)
    parent, source = predecessor(step, cfg, output)
    if parent:
        failure = first_invalid(parent, predecessor(parent, cfg, source)[1], cfg, source)
        if failure:
            return failure
    if not output.exists():
        return step, f"missing checkpoint {output}"
    if not manifest_path(output).exists():
        return step, f"missing validation manifest for {output}"
    try:
        saved = json.loads(manifest_path(output).read_text())
        if saved['output_sha256'] != digest(output):
            return step, f"checkpoint checksum changed: {output}; it may be incomplete or overwritten"
        if saved['signature'] != signature(step, input_path, cfg):
            return step, f"inputs, parameters or code changed for {output}"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return step, f"cannot verify {output}: {exc}"
    return None


def require_predecessor(step, cfg, output):
    parent, source = predecessor(step,cfg,output)
    failure = first_invalid(parent, predecessor(parent,cfg,source)[1], cfg, source) if parent else None
    if failure:
        failed_step, reason = failure
        raise ValueError(f"Cannot run {step}: {reason}. Rebuild from '{failed_step}' through '{parent}'.")


def atomic_write_h5ad(adata, output):
    """Finish HDF5 serialization before replacing the previous checkpoint."""
    import os
    import tempfile
    output = Path(output)
    fd, name = tempfile.mkstemp(prefix=f'.{output.name}.', suffix='.tmp.h5ad', dir=output.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        adata.write_h5ad(temporary)
        with temporary.open('rb') as stream:
            os.fsync(stream.fileno())
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)

def record(step, input_path, cfg, output):
    saved={'signature':signature(step,input_path,cfg),'output_sha256':digest(output)}
    destination=manifest_path(output); temporary=Path(str(destination)+'.tmp')
    temporary.write_text(json.dumps(saved,indent=2)); temporary.replace(destination)
