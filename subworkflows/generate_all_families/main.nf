#!/usr/bin/env nextflow

include { CREATE_CLUSTERS_PKL } from "${params.moduleDir}/family/main.nf"
include { REFINE_FAMILIES     } from "${params.moduleDir}/family/main.nf"

workflow GENERATE_FAMILIES {
    take:
    clusters_tsv
    mgnifams_fasta
    starting_num_sequences

    main:
    res              = CREATE_CLUSTERS_PKL(clusters_tsv)
    refined_families = REFINE_FAMILIES(res.pkl, res.refined_families, mgnifams_fasta, \
        res.discarded_clusters, res.converged_families, res.family_metadata, \
        starting_num_sequences, params.minimum_members, params.iteration)

    emit:
    tsv          = refined_families.tsv
    seed_msa_sto = refined_families.seed_msa_sto
    msa_sto      = refined_families.msa_sto
    hmm          = refined_families.hmm
    domtblout    = refined_families.domtblout
}
