#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.config   = "config/schema.yaml"
params.outdir   = "results"
params.modality = "scrna"
params.help     = false

if (params.help) {
    log.info """
    OmicSage Pipeline
    Usage:
        nextflow run pipeline/main.nf --modality scrna --outdir results/
    Modalities: scrna | scatac | spatial | multiome
    """.stripIndent()
    System.exit(0)
}

include { SCRNA_WORKFLOW    } from './workflows/scrna.nf'
include { SCATAC_WORKFLOW   } from './workflows/scatac.nf'
include { SPATIAL_WORKFLOW  } from './workflows/spatial.nf'
include { MULTIOME_WORKFLOW } from './workflows/integration.nf'

log.info """
OmicSage v0.1.0-dev
  Modality : ${params.modality}
  Config   : ${params.config}
  Output   : ${params.outdir}
""".stripIndent()

workflow {
    switch (params.modality) {
        case 'scrna':    SCRNA_WORKFLOW();    break
        case 'scatac':   SCATAC_WORKFLOW();   break
        case 'spatial':  SPATIAL_WORKFLOW();  break
        case 'multiome': MULTIOME_WORKFLOW(); break
        default: error "Unknown modality: '${params.modality}'"
    }
}

workflow.onComplete {
    log.info "Pipeline ${workflow.success ? 'COMPLETE ✓' : 'FAILED ✗'} — ${workflow.duration}"
}
