process CITE_HARMONY_ADT {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path reduce_checkpoint      // cite_03_reduced_adt.h5ad

    output:
    path "cite_04_harmony_adt.h5ad", emit: checkpoint

    script:
    """
    python /app/run_cite_pipeline.py \\
        --config ${config} \\
        --step harmony_adt
    """
}
