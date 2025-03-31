include { BUILD_PYFASTX_INDEX        } from "../../../modules/local/build_pyfastx_index/main"
include { GENERATE_FAMILIES          } from "../../../modules/local/generate_families/main"
include { REMOVE_REDUNDANCY          } from "../../../subworkflows/local/remove_redundancy"
include { PRESENT_DISCARDED_FAMILIES } from "../../../modules/local/present_discarded_families/main"
include { PRESENT_FAMILY_METADATA    } from "../../../modules/local/present_family_metadata/main"

workflow GENERATE_NONREDUNDANT_FAMILIES {
    take:
    cluster_chunks
    mgnifams_fa

    main:
    ch_versions = Channel.empty()

    ch_pyfastx = BUILD_PYFASTX_INDEX( mgnifams_fa )
    ch_versions = ch_versions.mix( BUILD_PYFASTX_INDEX.out.versions )

    ch_families = GENERATE_FAMILIES( cluster_chunks, mgnifams_fa.first(), ch_pyfastx.index.first() )
    ch_versions = ch_versions.mix( GENERATE_FAMILIES.out.versions )

    ch_hmm = ch_families.hmm
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"pre_redundant"], file ] }

    ch_reps_fasta = ch_families.fasta
        .map { meta, files -> files }
        .collect()
        .map{ file -> [[id: 'reps_fasta'], file] }

    ch_metadata = ch_families.metadata
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"metadata"], file ] }

    ch_seed_msa_sto = ch_families.seed_msa_sto
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"seed_msa_sto"], file ] }

    ch_msa_sto = ch_families.msa_sto
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"msa_sto"], file ] }

    ch_rf = ch_families.rf
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"rf"], file ] }

    ch_domtblout = ch_families.domtblout
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"domtblout"], file ] }

    ch_tsv = ch_families.tsv
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"tsv"], file ] }

    ch_discarded = ch_families.discarded
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"discarded"], file ] }

    ch_successful = ch_families.successful
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"successful"], file ] }

    ch_converged = ch_families.converged
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"converged"], file ] }

    ch_logs = ch_families.logs
        .map { meta, files -> files }
        .collect()
        .map { file -> [ [id:"logs"], file ] }

    REMOVE_REDUNDANCY( ch_hmm, ch_reps_fasta, ch_metadata, \
        ch_seed_msa_sto, ch_msa_sto, ch_rf, ch_domtblout, ch_tsv, \
        ch_discarded, ch_successful, ch_converged, ch_logs )
    ch_versions = ch_versions.mix( REMOVE_REDUNDANCY.out.versions )

    PRESENT_DISCARDED_FAMILIES( REMOVE_REDUNDANCY.out.discarded )
    ch_versions = ch_versions.mix( PRESENT_DISCARDED_FAMILIES.out.versions )

    PRESENT_FAMILY_METADATA( REMOVE_REDUNDANCY.out.metadata )
    ch_versions = ch_versions.mix( PRESENT_FAMILY_METADATA.out.versions )

    emit:
    versions      = ch_versions
    seed_msa_sto  = REMOVE_REDUNDANCY.out.seed_msa_sto
    msa_sto       = REMOVE_REDUNDANCY.out.msa_sto
    hmm           = REMOVE_REDUNDANCY.out.hmm
    rf            = REMOVE_REDUNDANCY.out.rf
    tsv           = REMOVE_REDUNDANCY.out.tsv
    converged     = REMOVE_REDUNDANCY.out.converged
    metadata      = REMOVE_REDUNDANCY.out.metadata
    family_reps   = REMOVE_REDUNDANCY.out.family_reps
    discarded_mqc = PRESENT_DISCARDED_FAMILIES.out.mqc
    metadata_mqc  = PRESENT_FAMILY_METADATA.out.mqc
}
