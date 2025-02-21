#!/usr/bin/env nextflow
// TODO remove
include { CHUNK_CLUSTERS           } from "../../../modules/local/chunk_clusters/main"
include { REFINE_FAMILIES_PARALLEL } from "../../../modules/local/refine_families_parallel/main"

workflow GENERATE_FAMILIES_PARALLEL {
    take:
    clusters_tsv
    checked_clusters
    mgnifams_fasta

    main:
    ch_versions = Channel.empty()

    cluster_chunks = CHUNK_CLUSTERS( clusters_tsv, checked_clusters, params.minimum_members, params.num_cluster_chunks ).fasta_chunks
    ch_versions = ch_versions.mix( CHUNK_CLUSTERS.out.versions )

    ch_cluster = cluster_chunks
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file_path.getSimpleName().split('_')[-1]], file_path ]
        }

    refined_families = REFINE_FAMILIES_PARALLEL( ch_cluster, mgnifams_fasta.first() )
    ch_versions = ch_versions.mix( REFINE_FAMILIES_PARALLEL.out.versions )

    emit:
    versions     = ch_versions
    seed_msa_sto = refined_families.seed_msa_sto
    msa_sto      = refined_families.msa_sto
    hmm          = refined_families.hmm
    rf           = refined_families.rf
    domtblout    = refined_families.domtblout
    tsv          = refined_families.tsv
    discarded    = refined_families.discarded
    successful   = refined_families.successful
    converged    = refined_families.converged
    metadata     = refined_families.metadata
    logs         = refined_families.logs
}
