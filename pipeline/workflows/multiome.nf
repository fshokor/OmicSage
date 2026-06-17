#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MULTIOME_ATAC_QC       } from '../modules/multiome/atac_qc.nf'
include { MULTIOME_ATAC_REDUCE   } from '../modules/multiome/atac_reduce.nf'
include { MULTIOME_ATAC_ANNOTATE } from '../modules/multiome/atac_annotate.nf'
include { MULTIOME_INTEGRATION   } from '../modules/multiome/multiome_integration.nf'
include { MULTIOME_DEG           } from '../modules/multiome/multiome_deg.nf'
include { MULTIOME_GRN           } from '../modules/multiome/multiome_grn.nf'

workflow MULTIOME_WORKFLOW {
    take:
    config_ch

    main:
    MULTIOME_ATAC_QC(config_ch)
    MULTIOME_ATAC_REDUCE(config_ch, MULTIOME_ATAC_QC.out.checkpoint)
    MULTIOME_ATAC_ANNOTATE(config_ch, MULTIOME_ATAC_REDUCE.out.checkpoint)
    MULTIOME_INTEGRATION(config_ch, MULTIOME_ATAC_ANNOTATE.out.checkpoint)
    MULTIOME_DEG(config_ch, MULTIOME_INTEGRATION.out.checkpoint)
    MULTIOME_GRN(config_ch, MULTIOME_DEG.out.checkpoint)

    emit:
    grn = MULTIOME_GRN.out.checkpoint
}
