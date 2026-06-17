process CITE_NORMALIZE_ADT {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config

    output:
    path "cite_01_normalized_adt.h5ad", emit: checkpoint

    script:
    """
    python /app/run_cite_pipeline.py \\
        --config ${config} \\
        --step normalize_adt
    """
}
