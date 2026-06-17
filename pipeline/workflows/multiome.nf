#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── OmicSage — workflows/multiome.nf ──────────────────────────────────────────
// Multiome workflow: chains all 6 pipeline steps.
// Checkpoint files: multiome_01_*.h5ad → multiome_06_*.h5mu
//
// Dependency note: atac_annotate and multiome_integration both also read
// rna_input from the config (paths.rna_input). Nextflow passes the config
// and the Python runner resolves rna_input internally.

include { MULTIOME_ATAC_QC            } from '../modules/multiome/atac_qc.nf'
include { MULTIOME_ATAC_REDUCE        } from '../modules/multiome/atac_reduce.nf'
include { MULTIOME_ATAC_ANNOTATE      } from '../modules/multiome/atac_annotate.nf'
include { MULTIOME_INTEGRATION        } from '../modules/multiome/multiome_integration.nf'
include { MULTIOME_DEG                } from '../modules/multiome/multiome_deg.nf'
include { MULTIOME_GRN                } from '../modules/multiome/multiome_grn.nf'

workflow MULTIOME_WORKFLOW {

    take:
    config_ch   // channel: path to run config YAML

    main:

    MULTIOME_ATAC_QC(config_ch)

    MULTIOME_ATAC_REDUCE(
        config_ch,
        MULTIOME_ATAC_QC.out.checkpoint
    )

    // atac_annotate also needs rna_input — resolved from config by the runner
    MULTIOME_ATAC_ANNOTATE(
        config_ch,
        MULTIOME_ATAC_REDUCE.out.checkpoint
    )

    // multiome_integration also needs rna_input — resolved from config
    MULTIOME_INTEGRATION(
        config_ch,
        MULTIOME_ATAC_ANNOTATE.out.checkpoint
    )

    MULTIOME_DEG(
        config_ch,
        MULTIOME_INTEGRATION.out.checkpoint
    )

    MULTIOME_GRN(
        config_ch,
        MULTIOME_DEG.out.checkpoint
    )

    emit:
    grn = MULTIOME_GRN.out.checkpoint
}
