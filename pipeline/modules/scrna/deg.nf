process SCRNA_DEG {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path annotate_checkpoint    // 05_annotated.h5ad

    output:
    path "06_deg.h5ad", emit: checkpoint

    script:
    """
    python /app/run_scrna_pipeline.py \\
        --config ${config} \\
        --step deg
    """
}
