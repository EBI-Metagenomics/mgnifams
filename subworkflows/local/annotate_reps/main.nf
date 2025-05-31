#!/usr/bin/env nextflow

include { S4PRED_RUNMODEL                } from '../../../modules/nf-core/s4pred/runmodel/main' 
include { PARSE_S4PRED_TO_FEATURE_VIEWER } from '../../../modules/local/parse_s4pred_to_feature_viewer/main' 

workflow ANNOTATE_REPS {
    take:
    fasta
    
    main:
    ch_versions = Channel.empty()

    S4PRED_RUNMODEL( fasta )
    ch_versions = ch_versions.mix( S4PRED_RUNMODEL.out.versions )

    PARSE_S4PRED_TO_FEATURE_VIEWER( S4PRED_RUNMODEL.out.preds )
    ch_versions = ch_versions.mix( PARSE_S4PRED_TO_FEATURE_VIEWER.out.versions )

    emit:
    versions        = ch_versions
    s4preds         = S4PRED_RUNMODEL.out.preds
    s4pred_features = PARSE_S4PRED_TO_FEATURE_VIEWER.out.features
}
