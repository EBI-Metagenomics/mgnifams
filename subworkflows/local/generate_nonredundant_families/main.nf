include { HMMER_ESLSFETCHINDEX                  } from '../../../modules/nf-core/hmmer/eslsfetchindex/main'
include { MGNIFAM_GENERATEFAMILIES              } from '../../../modules/nf-core/mgnifam/generatefamilies/main'
include { REMOVE_REDUNDANCY                     } from '../../../subworkflows/local/remove_redundancy'
include { HHSUITE_REFORMAT as REFORMAT_SEED_MSA } from '../../../modules/nf-core/hhsuite/reformat/main'
include { HHSUITE_REFORMAT as REFORMAT_FULL_MSA } from '../../../modules/nf-core/hhsuite/reformat/main'
include { PRESENT_DISCARDED_FAMILIES            } from '../../../modules/local/present_discarded_families/main'
include { PRESENT_FAMILY_METADATA               } from '../../../modules/local/present_family_metadata/main'

workflow GENERATE_NONREDUNDANT_FAMILIES {
    take:
    cluster_chunks // channel: [ meta, clusters_chunk.tsv ]; the mgnifams_* generation params reach mgnifam via ext.args
    mgnifams_fa    // channel: [ meta, mgnifams_input.fa ] (plain: the SSI index needs an uncompressed file)
    outdir
    redundant_length_threshold
    redundant_score_threshold
    similarity_score_threshold
    starting_id

    main:
    ch_versions = channel.empty()

    // One index of the input FASTA, shared by every chunk
    HMMER_ESLSFETCHINDEX( mgnifams_fa )
    ch_generate_input = cluster_chunks
        .combine( mgnifams_fa.join( HMMER_ESLSFETCHINDEX.out.ssi ).map { _meta, fa, ssi -> [ fa, ssi ] } )
    ch_families = MGNIFAM_GENERATEFAMILIES( ch_generate_input )

    ch_hmm         = collectedAs( ch_families.hmm       , 'pre_redundant' )
    ch_reps_fasta  = collectedAs( ch_families.reps_fasta, 'reps_fasta' )
    ch_metadata    = collectedAs( ch_families.csv       , 'metadata' )
    ch_seed_msa_sto = collectedAs( ch_families.seed_msa , 'seed_msa_sto' )
    ch_full_msa_sto = collectedAs( ch_families.full_msa , 'full_msa_sto' )
    ch_rf          = collectedAs( ch_families.rf        , 'rf' )
    ch_tsv         = collectedAs( ch_families.tsv       , 'tsv' )
    ch_discarded   = collectedAs( ch_families.discarded , 'discarded' )
    ch_successful  = collectedAs( ch_families.successful, 'successful' )
    ch_converged   = collectedAs( ch_families.converged , 'converged' )
    ch_logs        = collectedAs( ch_families.log       , 'logs' )

    REMOVE_REDUNDANCY( ch_reps_fasta, outdir, ch_hmm, ch_metadata, \
        ch_seed_msa_sto, ch_full_msa_sto, ch_rf, ch_tsv, \
        redundant_length_threshold, redundant_score_threshold, \
        similarity_score_threshold, ch_discarded, ch_successful, \
        ch_converged, ch_logs, starting_id )
    ch_versions = ch_versions.mix( REMOVE_REDUNDANCY.out.versions )

    ch_input_for_reformat_seed = REMOVE_REDUNDANCY.out.seed_msa_sto
        .transpose()
        .map { meta, file ->
            [[id: meta.id, chunk: file.getSimpleName()], file]
        }
    REFORMAT_SEED_MSA(ch_input_for_reformat_seed, "sto", "fas")

    ch_input_for_reformat_full = REMOVE_REDUNDANCY.out.full_msa_sto
        .transpose()
        .map { meta, file ->
            [[id: meta.id, chunk: file.getSimpleName()], file]
        }
    REFORMAT_FULL_MSA(ch_input_for_reformat_full, "sto", "fas")

    PRESENT_DISCARDED_FAMILIES( REMOVE_REDUNDANCY.out.discarded )
    ch_versions = ch_versions.mix( PRESENT_DISCARDED_FAMILIES.out.versions )

    PRESENT_FAMILY_METADATA( REMOVE_REDUNDANCY.out.metadata )
    ch_versions = ch_versions.mix( PRESENT_FAMILY_METADATA.out.versions )

    emit:
    versions         = ch_versions
    seed_msa_sto     = REMOVE_REDUNDANCY.out.seed_msa_sto
    full_msa_sto     = REMOVE_REDUNDANCY.out.full_msa_sto
    hmm              = REMOVE_REDUNDANCY.out.hmm
    rf               = REMOVE_REDUNDANCY.out.rf
    tsv              = REMOVE_REDUNDANCY.out.tsv
    converged        = REMOVE_REDUNDANCY.out.converged
    metadata         = REMOVE_REDUNDANCY.out.metadata
    family_reps      = REMOVE_REDUNDANCY.out.family_reps
    family_ids_fasta = REMOVE_REDUNDANCY.out.family_ids_fasta
    seed_msa         = REFORMAT_SEED_MSA.out.msa
    full_msa         = REFORMAT_FULL_MSA.out.msa
    discarded_mqc    = PRESENT_DISCARDED_FAMILIES.out.mqc
    metadata_mqc     = PRESENT_FAMILY_METADATA.out.mqc
    similarity_mqc   = REMOVE_REDUNDANCY.out.similarity_mqc
}

// All per-chunk files of one MGNIFAM_GENERATEFAMILIES output, as one [ [id: id], files ] item.
// Sorted by name: collect() keeps task-completion order, which would change the task hash and break -resume.
def collectedAs(ch, id) {
    return ch.map { _meta, files -> files }.collect(sort: { a, b -> a.name <=> b.name }).map { files -> [ [id: id], files ] }
}
