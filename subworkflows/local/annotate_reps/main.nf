#!/usr/bin/env nextflow

include { S4PRED_RUNMODEL } from '../../../modules/nf-core/s4pred/runmodel/main' 

workflow ANNOTATE_REPS {
    take:
    fasta
    
    main:
    ch_versions = Channel.empty()

    S4PRED_RUNMODEL( fasta )
    ch_versions = ch_versions.mix( S4PRED_RUNMODEL.out.versions )

    emit:
    versions = ch_versions
    s4preds  = S4PRED_RUNMODEL.out.preds
}
