/*
    SEQUENCE CLUSTERING
*/

include { EXTRACT_UNANNOTATED_FASTA } from "../../../subworkflows/local/extract_unannotated_fasta"
include { EXECUTE_CLUSTERING        } from "../../../subworkflows/local/execute_clustering"
include { PARSE_CLUSTER_STATS       } from "../../../modules/local/parse_cluster_stats/main"

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

    PARSE_CLUSTER_STATS( EXECUTE_CLUSTERING.out.clusters_tsv ) // multiqc output for cluster size distribution
    // TODO FILTER_CLUSTER_REP_NAMES( EXECUTE_CLUSTERING.out.clusters_tsv ) // awk probably
    // TODO splitText, cluster_size
    // TODO CHUNK_CLUSTERS( EXECUTE_CLUSTERING.out.clusters_tsv, ch_from_split )

    emit:
    versions          = ch_versions
    mgnifams_input_fa = ch_mgnifams_input_fa
    clusters_tsv      = EXECUTE_CLUSTERING.out.clusters_tsv
}
