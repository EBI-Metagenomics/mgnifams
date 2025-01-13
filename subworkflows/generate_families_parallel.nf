#!/usr/bin/env nextflow

include { CHUNK_CLUSTERS           } from "../modules/local/family/chunk_clusters.nf"
include { REFINE_FAMILIES_PARALLEL } from "../modules/local/family/refine_families_parallel.nf"

workflow GENERATE_FAMILIES_PARALLEL {
    take:
    clusters_tsv
    checked_clusters
    mgnifams_fasta

    main:
    clusters_chunks  = CHUNK_CLUSTERS(clusters_tsv, checked_clusters, params.minimum_members, params.num_cluster_chunks)
    
    //TODO
    CHUNK_CLUSTERS.out.view()
    // refined_families = REFINE_FAMILIES_PARALLEL(clusters_chunks.flatten(), mgnifams_fasta.first())
    
    // emit:
    // seed_msa_sto = refined_families.seed_msa_sto
    // msa_sto      = refined_families.msa_sto
    // hmm          = refined_families.hmm
    // rf           = refined_families.rf
    // domtblout    = refined_families.domtblout
    // tsv          = refined_families.tsv
    // discarded    = refined_families.discarded
    // successful   = refined_families.successful
    // converged    = refined_families.converged
    // metadata     = refined_families.metadata
    // logs         = refined_families.logs
}
