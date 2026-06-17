process MULTIOME_DEG {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path integration_checkpoint // multiome_04_integration.h5mu

    output:
    path "multiome_05_deg.h5mu", emit: checkpoint

    script:
    """
    python /app/run_multiome_pipeline.py \\
        --config ${config} \\
        --step multiome_deg
    """
}
