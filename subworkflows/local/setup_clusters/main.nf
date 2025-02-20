/*
    SEQUENCE CLUSTERING
*/

include { EXTRACT_UNANNOTATED_FASTA } from "../../../subworkflows/local/extract_unannotated_fasta"
include { EXECUTE_CLUSTERING        } from "../../../subworkflows/local/execute_clustering"

workflow SETUP_CLUSTERS {
    take:
    input

    main:
    ch_versions = Channel.empty()

    if (!params.fasta_input_mode) {
        ch_mgnifams_input_fa = EXTRACT_UNANNOTATED_FASTA( input ).fasta
        ch_versions = ch_versions.mix( EXTRACT_UNANNOTATED_FASTA.out.versions )
    } else {
        ch_mgnifams_input_fa = channel.fromPath(input)
    }
    EXECUTE_CLUSTERING( ch_mgnifams_input_fa )
    ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

    emit:
    versions          = ch_versions
    mgnifams_input_fa = ch_mgnifams_input_fa
    clusters_tsv      = EXECUTE_CLUSTERING.out.clusters_tsv
    num_sequences     = EXECUTE_CLUSTERING.out.num_sequences
}
