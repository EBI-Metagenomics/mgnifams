include { S4PRED_RUNMODEL                } from '../../../modules/nf-core/s4pred/runmodel/main'
include { PARSE_S4PRED_TO_FEATURE_VIEWER } from '../../../modules/local/parse_s4pred_to_feature_viewer/main'
include { DEEPTMHMM_PREDICT              } from '../../../modules/local/deeptmhmm/predict/main'
include { PARSE_TM_TO_FEATURE_VIEWER     } from '../../../modules/local/parse_tm_to_feature_viewer/main'
include { HMMER_HMMSEARCH                } from '../../../modules/nf-core/hmmer/hmmsearch/main'

workflow ANNOTATE_REPS {
    take:
    fasta
    skip_deeptmhmm
    deeptmhmm_path
    pfam_path
    funfams_path
    
    main:
    ch_versions       = Channel.empty()
    ch_tm_composition = Channel.of([ [ id: 'reps_fasta' ], [] ])

    S4PRED_RUNMODEL( fasta )
    ch_versions = ch_versions.mix( S4PRED_RUNMODEL.out.versions )

    PARSE_S4PRED_TO_FEATURE_VIEWER( S4PRED_RUNMODEL.out.preds )
    ch_versions = ch_versions.mix( PARSE_S4PRED_TO_FEATURE_VIEWER.out.versions )

    if (!skip_deeptmhmm && !workflow.profile.contains("conda")) {
        DEEPTMHMM_PREDICT( fasta, deeptmhmm_path )
        ch_versions = ch_versions.mix( DEEPTMHMM_PREDICT.out.versions )

        ch_tm_composition = PARSE_TM_TO_FEATURE_VIEWER( DEEPTMHMM_PREDICT.out.line3 ).composition
        ch_versions = ch_versions.mix( PARSE_TM_TO_FEATURE_VIEWER.out.versions )
    }

    // TODO ch_pfam

    ch_funfams = Channel.of([ [ id: 'reps_fasta' ], file(funfams_path, checkIfExists: true) ])
    ch_input_for_hmmsearch = ch_funfams
        .combine(fasta, by: 0)
        .map { meta, model, seqs -> [meta, model, seqs, false, false, true] }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    emit:
    versions        = ch_versions
    s4preds         = S4PRED_RUNMODEL.out.preds
    s4pred_features = PARSE_S4PRED_TO_FEATURE_VIEWER.out.features
    composition     = PARSE_S4PRED_TO_FEATURE_VIEWER.out.composition
    tm_composition  = ch_tm_composition
    funfam_domains  = HMMER_HMMSEARCH.out.domain_summary
}
