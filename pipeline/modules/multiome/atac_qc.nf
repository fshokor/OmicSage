process MULTIOME_ATAC_QC {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config

    output:
    path "multiome_01_qc_atac.h5ad", emit: checkpoint

    script:
    """
    python /app/run_multiome_pipeline.py \\
        --config ${config} \\
        --step atac_qc
    """
}
