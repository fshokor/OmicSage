process SCRNA_REDUCE {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path normalize_checkpoint   // 02_normalized.h5ad

    output:
    path "03_reduced.h5ad", emit: checkpoint

    script:
    """
    python /app/run_scrna_pipeline.py \\
        --config ${config} \\
        --step reduce
    """
}
