process SPATIAL_IMPUTE {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path cluster_checkpoint     // 04_clustered.h5ad (predecessor = cluster, not deconvolve)

    output:
    path "07_imputed.h5ad", emit: checkpoint

    script:
    """
    python /app/run_spatial_pipeline.py \\
        --config ${config} \\
        --step impute
    """
}
