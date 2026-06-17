#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── OmicSage — workflows/cite.nf ──────────────────────────────────────────────
// CITE-seq workflow: chains all 10 pipeline steps.
// Checkpoint files: cite_01_*.h5ad → cite_10_*.h5mu
//
// Dependency note: protein_rna_corr and epitope_characterisation both branch
// from integration (cite_06), not from gsea_cite, matching STEP_PREDECESSOR
// in run_cite_pipeline.py.

include { CITE_NORMALIZE_ADT          } from '../modules/cite/normalize_adt.nf'
include { CITE_DOUBLETS               } from '../modules/cite/doublets.nf'
include { CITE_REDUCE_ADT             } from '../modules/cite/reduce_adt.nf'
include { CITE_HARMONY_ADT            } from '../modules/cite/harmony_adt.nf'
include { CITE_ANNOTATE_ADT           } from '../modules/cite/annotate_adt.nf'
include { CITE_INTEGRATION            } from '../modules/cite/integration.nf'
include { CITE_DEG                    } from '../modules/cite/deg_cite.nf'
include { CITE_GSEA                   } from '../modules/cite/gsea_cite.nf'
include { CITE_PROTEIN_RNA_CORR       } from '../modules/cite/protein_rna_corr.nf'
include { CITE_EPITOPE_CHARACTERISATION } from '../modules/cite/epitope_characterisation.nf'

workflow CITE_WORKFLOW {

    take:
    config_ch   // channel: path to run config YAML

    main:

    // ── linear chain: normalize → doublets → reduce → harmony → annotate ─────
    CITE_NORMALIZE_ADT(config_ch)

    CITE_DOUBLETS(
        config_ch,
        CITE_NORMALIZE_ADT.out.checkpoint
    )

    CITE_REDUCE_ADT(
        config_ch,
        CITE_DOUBLETS.out.checkpoint
    )

    CITE_HARMONY_ADT(
        config_ch,
        CITE_REDUCE_ADT.out.checkpoint
    )

    CITE_ANNOTATE_ADT(
        config_ch,
        CITE_HARMONY_ADT.out.checkpoint
    )

    // integration reads annotate_adt output + rna_input (resolved from config)
    CITE_INTEGRATION(
        config_ch,
        CITE_ANNOTATE_ADT.out.checkpoint
    )

    // ── linear chain post-integration ────────────────────────────────────────
    CITE_DEG(
        config_ch,
        CITE_INTEGRATION.out.checkpoint
    )

    CITE_GSEA(
        config_ch,
        CITE_DEG.out.checkpoint
    )

    // ── branches off integration (not gsea) ──────────────────────────────────
    CITE_PROTEIN_RNA_CORR(
        config_ch,
        CITE_INTEGRATION.out.checkpoint
    )

    CITE_EPITOPE_CHARACTERISATION(
        config_ch,
        CITE_INTEGRATION.out.checkpoint
    )

    emit:
    gsea        = CITE_GSEA.out.checkpoint
    corr        = CITE_PROTEIN_RNA_CORR.out.checkpoint
    epitope     = CITE_EPITOPE_CHARACTERISATION.out.checkpoint
}
