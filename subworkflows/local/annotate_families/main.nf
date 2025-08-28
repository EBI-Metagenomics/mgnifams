include { ANNOTATE_REPS       } from '../../../subworkflows/local/annotate_reps'
include { ANNOTATE_MODELS     } from '../../../subworkflows/local/annotate_models'
include { ANNOTATE_STRUCTURES } from '../../../subworkflows/local/annotate_structures'

workflow ANNOTATE_FAMILIES {
    take:
    reps
    skip_deeptmhmm
    deeptmhmm_path
    pfam_path
    funfams_path
    seed_msa
    full_msa
    hh_mode
    hhdb_path
    pdb
    foldseek_db_path
    outdir

    main:
    ch_versions = Channel.empty()

    ANNOTATE_REPS( reps, skip_deeptmhmm, deeptmhmm_path, pfam_path, funfams_path )
    ch_versions = ch_versions.mix( ANNOTATE_REPS.out.versions )

    ANNOTATE_MODELS( seed_msa, hh_mode, hhdb_path )
    ch_versions = ch_versions.mix( ANNOTATE_MODELS.out.versions )

    ANNOTATE_STRUCTURES( pdb, foldseek_db_path, outdir )
    ch_versions = ch_versions.mix( ANNOTATE_STRUCTURES.out.versions )

    emit:
    versions           = ch_versions
    s4preds            = ANNOTATE_REPS.out.s4preds
    s4pred_features    = ANNOTATE_REPS.out.s4pred_features
    composition        = ANNOTATE_REPS.out.composition
    tm_composition     = ANNOTATE_REPS.out.tm_composition
    pfam_domains       = ANNOTATE_REPS.out.pfam_domains
    funfam_domains     = ANNOTATE_REPS.out.funfam_domains
    pfam_model_hits    = ANNOTATE_MODELS.out.pfam_hits
    foldseek_hits      = ANNOTATE_STRUCTURES.out.foldseek_hits
}
