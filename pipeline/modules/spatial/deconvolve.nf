process SPATIAL_DECONVOLVE {
    label 'process_python'
    tag  "deconvolve"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "05_deconvolved.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF2'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "spatial.deconvolve".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[deconvolve] disabled in config -- skipping')
    pathlib.Path('05_deconvolved.h5ad').touch()
    sys.exit(0)
PYEOF2
    /opt/conda/envs/omicsage/bin/python /app/run_spatial_pipeline.py \\
        --config /app/${config_path} \\
        --step deconvolve
    """
}
