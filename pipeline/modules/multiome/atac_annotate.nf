process MULTIOME_ATAC_ANNOTATE {
    label 'process_python'
    tag  "atac_annotate"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "multiome_03_annotate_atac.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python - << 'PYEOF2'
import yaml, pathlib, sys
cfg  = yaml.safe_load(open('/app/${config_path}'))
keys = "steps.atac_annotate".split('.')
node = cfg
for k in keys:
    node = node.get(k, {}) if isinstance(node, dict) else {}
enabled = node.get('enabled', True) if isinstance(node, dict) else True
if not enabled:
    print('[atac_annotate] disabled in config -- skipping')
    pathlib.Path('multiome_03_annotate_atac.h5ad').touch()
    sys.exit(0)
PYEOF2
    /opt/conda/envs/omicsage/bin/python /app/run_multiome_pipeline.py \\
        --config /app/${config_path} \\
        --step atac_annotate
    """
}
