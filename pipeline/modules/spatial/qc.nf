process SPATIAL_QC {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path ingest_checkpoint      // 01_ingested.h5ad

    output:
    path "02_qc.h5ad", emit: checkpoint

    script:
    """
    python /app/run_spatial_pipeline.py \\
        --config ${config} \\
        --step qc
    """
}
