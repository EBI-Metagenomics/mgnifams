include { POOL_PREREDUNDANT_FAMILIES_TSV } from '../../../modules/local/pool_preredundant_families_tsv/main'
include { MGNIFAM_COMPUTECLANS           } from '../../../modules/nf-core/mgnifam/computeclans/main'
include { SELECT_NONREDUNDANT_FAMILIES   } from '../../../modules/local/select_nonredundant_families/main'
include { POOL_NONREDUNDANT_FAMILIES     } from '../../../modules/local/pool_nonredundant_families/main'

workflow REMOVE_REDUNDANCY {
    take:
    mgnifams_fa // channel: [ meta, mgnifams_input.fa ]: lengths of whole-protein members
    reps_fasta
    hmm
    metadata
    seed_msa_sto
    full_msa_sto
    rf
    tsv
    redundant_score_threshold
    redundant_residue_threshold
    discarded
    successful
    converged
    logs
    starting_id

    main:
    ch_versions = channel.empty()

    // Members of every chunk under their pooling key <chunk>_<id>, so two chunks' local family 1 stay apart
    POOL_PREREDUNDANT_FAMILIES_TSV( tsv )
    ch_versions = ch_versions.mix( POOL_PREREDUNDANT_FAMILIES_TSV.out.versions )

    // Edges at residue Jaccard >= the threshold, plus containment-only edges the selector skips
    MGNIFAM_COMPUTECLANS(
        POOL_PREREDUNDANT_FAMILIES_TSV.out.tsv.combine( mgnifams_fa.map { _meta, fa -> fa } )
    )

    SELECT_NONREDUNDANT_FAMILIES( MGNIFAM_COMPUTECLANS.out.clan_edges, metadata, \
        redundant_score_threshold, redundant_residue_threshold )

    pooled_families = POOL_NONREDUNDANT_FAMILIES( seed_msa_sto, \
        full_msa_sto, hmm, rf, tsv, discarded, successful, \
        converged, metadata, reps_fasta, logs, SELECT_NONREDUNDANT_FAMILIES.out.txt, \
        starting_id )
    ch_versions = ch_versions.mix( POOL_NONREDUNDANT_FAMILIES.out.versions )

    emit:
    versions         = ch_versions
    seed_msa_sto     = pooled_families.seed_msa_sto
    full_msa_sto     = pooled_families.full_msa_sto
    hmm              = pooled_families.hmm
    rf               = pooled_families.rf
    tsv              = pooled_families.tsv
    discarded        = pooled_families.discarded
    successful       = pooled_families.successful
    converged        = pooled_families.converged
    metadata         = pooled_families.metadata
    family_reps      = pooled_families.family_reps
    family_ids_fasta = pooled_families.family_ids_fasta
    id_mapping       = pooled_families.id_mapping
}
