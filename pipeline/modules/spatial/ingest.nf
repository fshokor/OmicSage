process SPATIAL_INGEST {
    label 'process_python'
    tag  "ingest"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    output:
    val "01_ingested.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF2'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "spatial.ingest".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[ingest] disabled in config -- skipping')
    pathlib.Path('01_ingested.h5ad').touch()
    sys.exit(0)
PYEOF2
    /opt/conda/envs/omicsage/bin/python /app/run_spatial_pipeline.py \\
        --config /app/${config_path} \\
        --step ingest
    """
}
