include { CAT_CAT                    } from '../../../modules/nf-core/cat/cat'
include { HMMER_HMMSEARCH            } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { IDENTIFY_REDUNDANT_FAMS    } from '../../../modules/local/identify_redundant_fams/main'
include { POOL_NONREDUNDANT_FAMILIES } from '../../../modules/local/pool_nonredundant_families/main'

workflow REMOVE_REDUNDANCY {
    take:
    hmm
    reps_fasta
    metadata
    seed_msa_sto
    msa_sto
    rf
    domtblout
    tsv
    discarded
    successful
    converged
    logs

    main:
    ch_versions = Channel.empty()

    ch_reps_fasta = reps_fasta
        .map { meta, files -> files }
        .flatten()
        .collectFile(name: "pre_redundant_reps.fasta", storeDir: params.outdir + "/generate_families")
        .map { file -> [[id: 'pre_redundant'], file] }

    CAT_CAT( hmm )
    ch_versions = ch_versions.mix( CAT_CAT.out.versions )

    ch_input_for_hmmsearch = CAT_CAT.out.file_out
        .combine(ch_reps_fasta, by: 0)
        .map { meta, model, seqs -> [meta, model, seqs, false, false, true] }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    // TODO cleverer way to remove instead of hmmsearch results, e.g. Jaccard indices
    IDENTIFY_REDUNDANT_FAMS( HMMER_HMMSEARCH.out.domain_summary, metadata, params.redundant_length_threshold )
    ch_versions = ch_versions.mix( IDENTIFY_REDUNDANT_FAMS.out.versions )

    pooled_families = POOL_NONREDUNDANT_FAMILIES( seed_msa_sto, \
        msa_sto, hmm, rf, domtblout, tsv, \
        discarded, successful, converged, \
        metadata, reps_fasta, logs, IDENTIFY_REDUNDANT_FAMS.out.txt, params.starting_id )
    ch_versions = ch_versions.mix( POOL_NONREDUNDANT_FAMILIES.out.versions )

    emit:
    versions     = ch_versions
    seed_msa_sto = pooled_families.seed_msa_sto
    msa_sto      = pooled_families.msa_sto
    hmm          = pooled_families.hmm
    rf           = pooled_families.rf
    domtblout    = pooled_families.domtblout
    tsv          = pooled_families.tsv
    discarded    = pooled_families.discarded
    successful   = pooled_families.successful
    converged    = pooled_families.converged
    metadata     = pooled_families.metadata
    family_reps  = pooled_families.family_reps
    id_mapping   = pooled_families.id_mapping
}
