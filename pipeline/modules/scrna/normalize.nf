process SCRNA_NORMALIZE {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path qc_checkpoint      // 01_qc.h5ad — declares dependency, not passed to script

    output:
    path "02_normalized.h5ad", emit: checkpoint

    script:
    """
    python /app/run_scrna_pipeline.py \\
        --config ${config} \\
        --step normalize
    """
}
