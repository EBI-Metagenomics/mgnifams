#!/usr/bin/env nextflow

include { CREATE_CLUSTERS_PKL } from "${params.moduleDir}/family/main.nf"
include { REFINE_FAMILIES     } from "${params.moduleDir}/family/main.nf"

workflow GENERATE_FAMILIES {
    take:
    clusters_tsv
    refined_families_tsv
    mgnifams_fasta
    discarded_clusters
    converged_families
    family_sizes
    starting_num_sequences
    iteration

    main:
    clusters_pkl     = CREATE_CLUSTERS_PKL(clusters_tsv).pkl
    refined_families = REFINE_FAMILIES(clusters_pkl, refined_families_tsv, mgnifams_fasta, discarded_clusters, converged_families, family_sizes, \
        starting_num_sequences, params.minimum_members, iteration)

    emit:
    tsv          = refined_families.tsv
    seed_msa_sto = refined_families.seed_msa_sto
    msa_sto      = refined_families.msa_sto
    hmm          = refined_families.hmm
    domtblout    = refined_families.domtblout
}
