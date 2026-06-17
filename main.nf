#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCRNA_WORKFLOW    } from './pipeline/workflows/scrna.nf'
include { CITE_WORKFLOW     } from './pipeline/workflows/cite.nf'
include { MULTIOME_WORKFLOW } from './pipeline/workflows/multiome.nf'
include { SPATIAL_WORKFLOW  } from './pipeline/workflows/spatial.nf'

workflow {

    if (params.help) {
        log.info """
        ╔══════════════════════════════════════════════════════╗
        ║              OmicSage — Nextflow DSL2                ║
        ╚══════════════════════════════════════════════════════╝

        Usage:
          nextflow run main.nf [options]

        Required:
          --config    Path to run config YAML
          --modality  scrna | cite | multiome | spatial

        Optional:
          --outdir    Results directory (default: results)
          --help      Show this message

        Resume (skip completed steps):
          nextflow run main.nf --config <yaml> --modality <mod> -resume

        Examples:
          nextflow run main.nf --config config/runs/GSE166635.yaml --modality scrna -resume
          nextflow run main.nf --config config/runs/GSE194122_cite.yaml --modality cite -resume
          nextflow run main.nf --config config/runs/GSE194122_atac_multiome.yaml --modality multiome -resume
          nextflow run main.nf --config config/runs/kuppe_heart.yaml --modality spatial -resume
        """.stripIndent()
        exit 0
    }

    if (!params.config) {
        error "ERROR: --config is required. Example: --config config/runs/GSE166635.yaml"
    }

    // Pass config as a string value — NOT as a staged path.
    // Modules prepend /app/ to get the absolute path inside the Docker container.
    config_ch = Channel.value(params.config)

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
        CITE_WORKFLOW(config_ch)
    } else if (params.modality == 'multiome') {
        MULTIOME_WORKFLOW(config_ch)
    } else if (params.modality == 'spatial') {
        SPATIAL_WORKFLOW(config_ch)
    } else {
        error "Unknown modality '${params.modality}'. Valid: scrna, cite, multiome, spatial"
    }
}
