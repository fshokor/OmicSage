#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CITE_NORMALIZE_ADT            } from '../modules/cite/normalize_adt.nf'
include { CITE_DOUBLETS                 } from '../modules/cite/doublets.nf'
include { CITE_REDUCE_ADT               } from '../modules/cite/reduce_adt.nf'
include { CITE_HARMONY_ADT              } from '../modules/cite/harmony_adt.nf'
include { CITE_ANNOTATE_ADT             } from '../modules/cite/annotate_adt.nf'
include { CITE_INTEGRATION              } from '../modules/cite/integration.nf'
include { CITE_DEG                      } from '../modules/cite/deg_cite.nf'
include { CITE_GSEA                     } from '../modules/cite/gsea_cite.nf'
include { CITE_PROTEIN_RNA_CORR         } from '../modules/cite/protein_rna_corr.nf'
include { CITE_EPITOPE_CHARACTERISATION } from '../modules/cite/epitope_characterisation.nf'

workflow CITE_WORKFLOW {
    take:
    config_ch

    main:
    CITE_NORMALIZE_ADT(config_ch)
    CITE_DOUBLETS(config_ch, CITE_NORMALIZE_ADT.out.checkpoint)
    CITE_REDUCE_ADT(config_ch, CITE_DOUBLETS.out.checkpoint)
    CITE_HARMONY_ADT(config_ch, CITE_REDUCE_ADT.out.checkpoint)
    CITE_ANNOTATE_ADT(config_ch, CITE_HARMONY_ADT.out.checkpoint)
    CITE_INTEGRATION(config_ch, CITE_ANNOTATE_ADT.out.checkpoint)
    CITE_DEG(config_ch, CITE_INTEGRATION.out.checkpoint)
    CITE_GSEA(config_ch, CITE_DEG.out.checkpoint)
    CITE_PROTEIN_RNA_CORR(config_ch, CITE_INTEGRATION.out.checkpoint)
    CITE_EPITOPE_CHARACTERISATION(config_ch, CITE_INTEGRATION.out.checkpoint)

    emit:
    gsea    = CITE_GSEA.out.checkpoint
    corr    = CITE_PROTEIN_RNA_CORR.out.checkpoint
    epitope = CITE_EPITOPE_CHARACTERISATION.out.checkpoint
}
