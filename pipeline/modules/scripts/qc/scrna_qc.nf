nextflow.enable.dsl = 2

process QC_SCRNA {
    label 'process_python'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path anndata

    output:
    path "*.qc.h5ad",       emit: filtered
    path "*.qc_metrics.json", emit: metrics

    script:
    """
    echo 'QC_SCRNA stub — implement in Phase 1'
    """
}
