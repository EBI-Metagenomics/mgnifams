/*
    SEQUENCE CLUSTERING
*/

include { EXTRACT_UNANNOTATED_FASTA      } from "../../../subworkflows/local/extract_unannotated_fasta"
include { EXECUTE_CLUSTERING             } from "../../../subworkflows/local/execute_clustering"
include { CALCULATE_CLUSTER_DISTRIBUTION } from "../../../modules/local/calculate_cluster_distribution/main"
include { EXTRACT_UNIQUE_CLUSTER_REPS    } from "../../../modules/local/extract_unique_cluster_reps/main"
include { CHUNK_CLUSTERS                 } from "../../../modules/local/chunk_clusters/main"

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

    EXTRACT_UNIQUE_CLUSTER_REPS( EXECUTE_CLUSTERING.out.clusters_tsv, params.minimum_members )
    ch_versions = ch_versions.mix( EXTRACT_UNIQUE_CLUSTER_REPS.out.versions )

    ch_cluster_reps_chunks = EXTRACT_UNIQUE_CLUSTER_REPS.out.reps
        .splitText(file:true, by: params.clusters_chunk_size)
        .map { meta, file ->
            [[id: meta.id, chunk: file.getBaseName(1).split('\\.')[-1]], file]
        }

    CHUNK_CLUSTERS( ch_cluster_reps_chunks, EXECUTE_CLUSTERING.out.clusters_tsv.first() )
    ch_versions = ch_versions.mix( CHUNK_CLUSTERS.out.versions )

    emit:
    versions          = ch_versions
    mgnifams_input_fa = ch_mgnifams_input_fa
    cluster_distr_mqc = CALCULATE_CLUSTER_DISTRIBUTION.out.mqc
    cluster_chunks    = CHUNK_CLUSTERS.out.tsv
}
