process CITE_GSEA {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path deg_checkpoint         // cite_07_deg.h5mu

    output:
    path "cite_08_gsea.h5mu", emit: checkpoint

    script:
    """
    python /app/run_cite_pipeline.py \\
        --config ${config} \\
        --step gsea_cite
    """
}
