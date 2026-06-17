process MULTIOME_GRN {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path deg_checkpoint         // multiome_05_deg.h5mu

    output:
    path "multiome_06_grn.h5mu", emit: checkpoint

    script:
    """
    python /app/run_multiome_pipeline.py \\
        --config ${config} \\
        --step multiome_grn
    """
}
