#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── OmicSage — workflows/scrna.nf ─────────────────────────────────────────────
// scRNA-seq workflow: chains all 10 pipeline steps.
// Each process calls run_scrna_pipeline.py --step <name> inside the container.
// Step outputs (h5ad checkpoint files) flow through channels as dependencies.

include { SCRNA_QC            } from '../modules/scrna/qc.nf'
include { SCRNA_NORMALIZE     } from '../modules/scrna/normalize.nf'
include { SCRNA_REDUCE        } from '../modules/scrna/reduce.nf'
include { SCRNA_CLUSTER       } from '../modules/scrna/cluster.nf'
include { SCRNA_ANNOTATE      } from '../modules/scrna/annotate.nf'
include { SCRNA_DEG           } from '../modules/scrna/deg.nf'
include { SCRNA_GSEA          } from '../modules/scrna/gsea.nf'
include { SCRNA_HARMONY       } from '../modules/scrna/harmony.nf'
include { SCRNA_CLUSTER_HARMONY } from '../modules/scrna/cluster_harmony.nf'
include { SCRNA_PSEUDOBULK    } from '../modules/scrna/pseudobulk.nf'

workflow SCRNA_WORKFLOW {

    take:
    config_ch   // channel: path to run config YAML

    main:

    // ── linear chain ──────────────────────────────────────────────────────────
    // Each step receives the config + the output of its predecessor.
    // Nextflow uses these dependencies to:
    //   1. Enforce execution order
    //   2. Cache completed steps on -resume

    SCRNA_QC(config_ch)

    SCRNA_NORMALIZE(
        config_ch,
        SCRNA_QC.out.checkpoint
    )

    SCRNA_REDUCE(
        config_ch,
        SCRNA_NORMALIZE.out.checkpoint
    )

    SCRNA_CLUSTER(
        config_ch,
        SCRNA_REDUCE.out.checkpoint
    )

    SCRNA_ANNOTATE(
        config_ch,
        SCRNA_CLUSTER.out.checkpoint
    )

    SCRNA_DEG(
        config_ch,
        SCRNA_ANNOTATE.out.checkpoint
    )

    SCRNA_GSEA(
        config_ch,
        SCRNA_DEG.out.checkpoint
    )

    SCRNA_HARMONY(
        config_ch,
        SCRNA_GSEA.out.checkpoint
    )

    // cluster_harmony and pseudobulk both depend on their own predecessors
    // pseudobulk runs off annotated (step 5), not harmony — matching STEP_PREDECESSOR
    SCRNA_CLUSTER_HARMONY(
        config_ch,
        SCRNA_HARMONY.out.checkpoint
    )

    SCRNA_PSEUDOBULK(
        config_ch,
        SCRNA_ANNOTATE.out.checkpoint   // pseudobulk predecessor = annotate
    )

    emit:
    pseudobulk     = SCRNA_PSEUDOBULK.out.checkpoint
    cluster_harmony = SCRNA_CLUSTER_HARMONY.out.checkpoint
}
