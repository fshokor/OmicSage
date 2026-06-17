process SCRNA_HARMONY {
    label 'process_python'
    tag  "harmony"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "08_harmony.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "steps.harmony".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[harmony] disabled in config -- skipping')
    pathlib.Path('08_harmony.h5ad').touch()
    sys.exit(0)
PYEOF
    /opt/conda/envs/omicsage/bin/python /app/run_scrna_pipeline.py \\
        --config /app/${config_path} \\
        --step harmony
    """
}
