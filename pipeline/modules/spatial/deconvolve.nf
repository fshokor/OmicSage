process SPATIAL_DECONVOLVE {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path cluster_checkpoint     // 04_clustered.h5ad

    output:
    path "05_deconvolved.h5ad", emit: checkpoint

    script:
    """
    python /app/run_spatial_pipeline.py \\
        --config ${config} \\
        --step deconvolve
    """
}
