process CITE_PROTEIN_RNA_CORR {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path integration_checkpoint // cite_06_integration.h5mu (predecessor = integration)

    output:
    path "cite_09_corr.h5mu", emit: checkpoint

    script:
    """
    python /app/run_cite_pipeline.py \\
        --config ${config} \\
        --step protein_rna_corr
    """
}
