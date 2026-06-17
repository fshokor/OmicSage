process CITE_INTEGRATION {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path annotate_checkpoint    // cite_05_annotated_adt.h5ad

    output:
    path "cite_06_integration.h5mu", emit: checkpoint

    script:
    """
    python /app/run_cite_pipeline.py \\
        --config ${config} \\
        --step integration
    """
}
