#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SPATIAL_INGEST     } from '../modules/spatial/ingest.nf'
include { SPATIAL_QC         } from '../modules/spatial/qc.nf'
include { SPATIAL_REDUCE     } from '../modules/spatial/reduce.nf'
include { SPATIAL_CLUSTER    } from '../modules/spatial/cluster.nf'
include { SPATIAL_DECONVOLVE } from '../modules/spatial/deconvolve.nf'
include { SPATIAL_DOWNSTREAM } from '../modules/spatial/downstream.nf'
include { SPATIAL_IMPUTE     } from '../modules/spatial/impute.nf'

workflow SPATIAL_WORKFLOW {
    take:
    config_ch

    main:
    SPATIAL_INGEST(config_ch)
    SPATIAL_QC(config_ch, SPATIAL_INGEST.out.checkpoint)
    SPATIAL_REDUCE(config_ch, SPATIAL_QC.out.checkpoint)
    SPATIAL_CLUSTER(config_ch, SPATIAL_REDUCE.out.checkpoint)
    SPATIAL_DECONVOLVE(config_ch, SPATIAL_CLUSTER.out.checkpoint)
    SPATIAL_DOWNSTREAM(config_ch, SPATIAL_DECONVOLVE.out.checkpoint)
    SPATIAL_IMPUTE(config_ch, SPATIAL_CLUSTER.out.checkpoint)

    emit:
    downstream = SPATIAL_DOWNSTREAM.out.checkpoint
    impute     = SPATIAL_IMPUTE.out.checkpoint
}
