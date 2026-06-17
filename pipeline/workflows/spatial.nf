#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── OmicSage — workflows/spatial.nf ───────────────────────────────────────────
// Spatial transcriptomics workflow: chains all 7 pipeline steps.
// Checkpoint files: 01_ingested.h5ad → 07_imputed.h5ad
//
// Dependency notes (from STEP_PREDECESSOR in run_spatial_pipeline.py):
//   downstream prefers deconvolve checkpoint; falls back to cluster
//             (the Python runner handles the fallback — Nextflow declares
//              deconvolve as the dependency to enforce ordering)
//   impute    reads from cluster checkpoint (does NOT require deconvolve)

include { SPATIAL_INGEST      } from '../modules/spatial/ingest.nf'
include { SPATIAL_QC          } from '../modules/spatial/qc.nf'
include { SPATIAL_REDUCE      } from '../modules/spatial/reduce.nf'
include { SPATIAL_CLUSTER     } from '../modules/spatial/cluster.nf'
include { SPATIAL_DECONVOLVE  } from '../modules/spatial/deconvolve.nf'
include { SPATIAL_DOWNSTREAM  } from '../modules/spatial/downstream.nf'
include { SPATIAL_IMPUTE      } from '../modules/spatial/impute.nf'

workflow SPATIAL_WORKFLOW {

    take:
    config_ch   // channel: path to run config YAML

    main:

    SPATIAL_INGEST(config_ch)

    SPATIAL_QC(
        config_ch,
        SPATIAL_INGEST.out.checkpoint
    )

    SPATIAL_REDUCE(
        config_ch,
        SPATIAL_QC.out.checkpoint
    )

    SPATIAL_CLUSTER(
        config_ch,
        SPATIAL_REDUCE.out.checkpoint
    )

    SPATIAL_DECONVOLVE(
        config_ch,
        SPATIAL_CLUSTER.out.checkpoint
    )

    // downstream prefers deconvolve output; runner falls back to cluster
    SPATIAL_DOWNSTREAM(
        config_ch,
        SPATIAL_DECONVOLVE.out.checkpoint
    )

    // impute branches from cluster (independent of deconvolve)
    SPATIAL_IMPUTE(
        config_ch,
        SPATIAL_CLUSTER.out.checkpoint
    )

    emit:
    downstream = SPATIAL_DOWNSTREAM.out.checkpoint
    impute     = SPATIAL_IMPUTE.out.checkpoint
}
