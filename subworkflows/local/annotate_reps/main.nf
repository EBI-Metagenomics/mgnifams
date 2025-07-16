include { S4PRED_RUNMODEL                } from '../../../modules/nf-core/s4pred/runmodel/main'
include { PARSE_S4PRED_TO_FEATURE_VIEWER } from '../../../modules/local/parse_s4pred_to_feature_viewer/main'
include { HMMER_HMMSEARCH                } from '../../../modules/nf-core/hmmer/hmmsearch/main'

workflow ANNOTATE_REPS {
    take:
    fasta
    funfams_path
    
    main:
    ch_versions = Channel.empty()

    S4PRED_RUNMODEL( fasta )
    ch_versions = ch_versions.mix( S4PRED_RUNMODEL.out.versions )

    PARSE_S4PRED_TO_FEATURE_VIEWER( S4PRED_RUNMODEL.out.preds )
    ch_versions = ch_versions.mix( PARSE_S4PRED_TO_FEATURE_VIEWER.out.versions )

    ch_funfams = Channel.of([ [ id: 'reps_fasta' ], file(funfams_path, checkIfExists: true) ])
    ch_input_for_hmmsearch = ch_funfams
        .combine(fasta, by: 0)
        .map { meta, model, seqs -> [meta, model, seqs, false, false, true] }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    emit:
    versions           = ch_versions
    s4preds            = S4PRED_RUNMODEL.out.preds
    s4pred_features    = PARSE_S4PRED_TO_FEATURE_VIEWER.out.features
    s4pred_composition = PARSE_S4PRED_TO_FEATURE_VIEWER.out.composition
    funfams_domains    = HMMER_HMMSEARCH.out.domain_summary
}
