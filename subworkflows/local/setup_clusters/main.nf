/*
    SEQUENCE CLUSTERING
*/

include { EXTRACT_UNANNOTATED_FASTA      } from "../../../subworkflows/local/extract_unannotated_fasta"
include { EXECUTE_CLUSTERING             } from "../../../subworkflows/local/execute_clustering"
include { CALCULATE_CLUSTER_DISTRIBUTION } from "../../../modules/local/calculate_cluster_distribution/main"

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

    CALCULATE_CLUSTER_DISTRIBUTION( EXECUTE_CLUSTERING.out.clusters_tsv )
    ch_versions = ch_versions.mix( CALCULATE_CLUSTER_DISTRIBUTION.out.versions )

    // TODO FILTER_CLUSTER_REP_NAMES( EXECUTE_CLUSTERING.out.clusters_tsv ) // awk probably
    // TODO splitText, cluster_size
    // TODO CHUNK_CLUSTERS( EXECUTE_CLUSTERING.out.clusters_tsv, ch_from_split )

    emit:
    versions          = ch_versions
    mgnifams_input_fa = ch_mgnifams_input_fa
    cluster_distr_mqc = CALCULATE_CLUSTER_DISTRIBUTION.out.mqc
    // clusters_tsv      = EXECUTE_CLUSTERING.out.clusters_tsv // TODO remove probably
}
