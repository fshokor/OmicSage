process MULTIOME_ATAC_ANNOTATE {
    label 'process_python'
    tag  "${config.simpleName}"

    errorStrategy { task.exitStatus in [130, 137, 139] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    path config
    path reduce_checkpoint      // multiome_02_reduce_atac.h5ad
                                // also reads paths.rna_input from config

    output:
    path "multiome_03_annotate_atac.h5ad", emit: checkpoint

    script:
    """
    python /app/run_multiome_pipeline.py \\
        --config ${config} \\
        --step atac_annotate
    """
}
