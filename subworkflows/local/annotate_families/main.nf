include { ANNOTATE_REPS       } from '../../../subworkflows/local/annotate_reps'
include { ANNOTATE_MODELS     } from '../../../subworkflows/local/annotate_models'
include { ANNOTATE_STRUCTURES } from '../../../subworkflows/local/annotate_structures'

workflow ANNOTATE_FAMILIES {
    take:
    reps
    funfams_path
    seed_msa
    full_msa
    hh_mode
    hhdb_path
    // pdb
    // foldseek_db_path
    // outdir

    main:
    ch_versions = Channel.empty()

    ANNOTATE_REPS( reps, funfams_path )
    ch_versions = ch_versions.mix( ANNOTATE_REPS.out.versions )

    ANNOTATE_MODELS( seed_msa, hh_mode, hhdb_path )
    ch_versions = ch_versions.mix( ANNOTATE_MODELS.out.versions )

    // TODO continue
    // ANNOTATE_STRUCTURES( pdb )
    // ch_versions = ch_versions.mix( ANNOTATE_STRUCTURES.out.versions )


    // old, TODO remove
    // seed_msa
    //     .map { meta, files ->
    //         String filePath = files[0]
    //         int lastIndex = filePath.lastIndexOf('/')
    //         String seed_msa_dir = filePath.substring(0, lastIndex + 1)
    //         [ [id:"seed_msa"], file(seed_msa_dir) ]
    //     }
    //     .set { seed_msa_ch }
    
    // full_msa
    //     .map { meta, files ->
    //         String filePath = files[0]
    //         int lastIndex = filePath.lastIndexOf('/')
    //         String msa_dir = filePath.substring(0, lastIndex + 1)
    //         [ [id:"msa"], file(msa_dir) ]
    //     }
    //     .set { hmmalign_msa_ch }

    // ch_structures = PREDICT_STRUCTURES(hmmalign_msa_ch)
    // ch_versions = ch_versions.mix( PREDICT_STRUCTURES.out.versions )

    // ch_pdb    = ch_structures.pdb
    // ch_scores = ch_structures.scores
    // ch_cif    = ch_structures.cif


    emit:
    versions  = ch_versions
    pfam_hits = ANNOTATE_MODELS.out.pfam_hits
    // foldseek_hits = ANNOTATE_STRUCTURES.out.foldseek_hits
    // scores        = ch_scores
    // cif           = ch_cif
}
