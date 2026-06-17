process MULTIOME_GRN {
    label 'process_python'
    tag  "multiome_grn"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "multiome_06_grn.h5mu", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF2'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "steps.multiome_grn".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[multiome_grn] disabled in config -- skipping')
    pathlib.Path('multiome_06_grn.h5mu').touch()
    sys.exit(0)
PYEOF2
    /opt/conda/envs/omicsage/bin/python /app/run_multiome_pipeline.py \\
        --config /app/${config_path} \\
        --step multiome_grn
    """
}
