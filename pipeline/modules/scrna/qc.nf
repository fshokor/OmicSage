process SCRNA_QC {
    label 'process_python'
    tag  "${config.simpleName}"

    // Retry with more memory on failure
    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config

    output:
    path "01_qc.h5ad", emit: checkpoint

    script:
    """
    python /app/run_scrna_pipeline.py \\
        --config ${config} \\
        --step qc
    """
}
