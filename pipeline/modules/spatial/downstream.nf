process SPATIAL_DOWNSTREAM {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path deconvolve_checkpoint  // 05_deconvolved.h5ad (falls back to 04_clustered.h5ad in runner)

    output:
    path "06_downstream.h5ad", emit: checkpoint

    script:
    """
    python /app/run_spatial_pipeline.py \\
        --config ${config} \\
        --step downstream
    """
}
