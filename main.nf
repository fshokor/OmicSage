#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── OmicSage — main.nf ────────────────────────────────────────────────────────
// Entry point. Routes to the correct modality workflow.
// Usage:
//   nextflow run main.nf --config config/runs/GSE166635.yaml --modality scrna
//   nextflow run main.nf --config config/runs/GSE166635.yaml --modality scrna -resume

include { SCRNA_WORKFLOW } from './pipeline/workflows/scrna.nf'

// ── entry workflow ────────────────────────────────────────────────────────────
workflow {

    // ── help ─────────────────────────────────────────────────────────────────
    if (params.help) {
        log.info """
        ╔══════════════════════════════════════════════════════╗
        ║              OmicSage — Nextflow DSL2                ║
        ╚══════════════════════════════════════════════════════╝

        Usage:
          nextflow run main.nf [options]

        Required:
          --config    Path to run config YAML  (e.g. config/runs/GSE166635.yaml)
          --modality  Pipeline modality        (scrna | cite | multiome | spatial)

        Optional:
          --outdir    Results directory        (default: results)
          --help      Show this message

        Resume a previous run (skip completed steps):
          nextflow run main.nf --config config/runs/GSE166635.yaml --modality scrna -resume

        Profiles:
          -profile local        Local execution with Docker  (default)
          -profile slurm        SLURM HPC with Singularity
          -profile test         Small test run
        """.stripIndent()
        exit 0
    }

    // ── validate required params ──────────────────────────────────────────────
    if (!params.config) {
        error "ERROR: --config is required. Example: --config config/runs/GSE166635.yaml"
    }

    // Pass config file as a channel so each process receives it as a path
    config_ch = Channel.fromPath(params.config, checkIfExists: true)

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║              OmicSage — Nextflow DSL2                ║
    ╚══════════════════════════════════════════════════════╝
    modality  : ${params.modality}
    config    : ${params.config}
    outdir    : ${params.outdir}
    """.stripIndent()

    if (params.modality == 'scrna') {
        SCRNA_WORKFLOW(config_ch)
    } else if (params.modality == 'cite') {
        error "CITE-seq workflow not yet implemented. Coming in Phase 11."
    } else if (params.modality == 'multiome') {
        error "Multiome workflow not yet implemented. Coming in Phase 11."
    } else if (params.modality == 'spatial') {
        error "Spatial workflow not yet implemented. Coming in Phase 11."
    } else {
        error "Unknown modality '${params.modality}'. Valid: scrna, cite, multiome, spatial"
    }
}
