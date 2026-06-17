process SCRNA_NORMALIZE {
    label 'process_python'
    tag  "normalize"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "02_normalized.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "steps.normalize".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[normalize] disabled in config -- skipping')
    pathlib.Path('02_normalized.h5ad').touch()
    sys.exit(0)
PYEOF
    /opt/conda/envs/omicsage/bin/python /app/run_scrna_pipeline.py \\
        --config /app/${config_path} \\
        --step normalize
    """
}
