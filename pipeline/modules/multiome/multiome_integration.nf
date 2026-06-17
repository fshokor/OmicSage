process MULTIOME_INTEGRATION {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path annotate_checkpoint    // multiome_03_annotate_atac.h5ad
                                // also reads paths.rna_input from config

    output:
    path "multiome_04_integration.h5mu", emit: checkpoint

    script:
    """
    python /app/run_multiome_pipeline.py \\
        --config ${config} \\
        --step multiome_integration
    """
}
