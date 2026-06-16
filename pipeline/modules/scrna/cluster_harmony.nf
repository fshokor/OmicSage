process SCRNA_CLUSTER_HARMONY {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path harmony_checkpoint     // 08_harmony.h5ad

    output:
    path "09_harmony_clustered.h5ad", emit: checkpoint

    script:
    """
    python /app/run_scrna_pipeline.py \\
        --config ${config} \\
        --step cluster_harmony
    """
}
