process SCRNA_DEG {
    label 'process_python'
    tag  "deg"
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2
    input:
    val config_path
    val predecessor
    output:
    val "06_deg.h5ad", emit: checkpoint
    script:
    """
    /opt/conda/envs/omicsage/bin/python /app/run_scrna_pipeline.py \
        --config /app/${config_path} \
        --step deg
    """
}
