process CITE_REDUCE_ADT {
    label 'process_python'
    tag  "reduce_adt"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "cite_03_reduced_adt.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF2'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "steps.reduce_adt".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[reduce_adt] disabled in config -- skipping')
    pathlib.Path('cite_03_reduced_adt.h5ad').touch()
    sys.exit(0)
PYEOF2
    /opt/conda/envs/omicsage/bin/python /app/run_cite_pipeline.py \\
        --config /app/${config_path} \\
        --step reduce_adt
    """
}
