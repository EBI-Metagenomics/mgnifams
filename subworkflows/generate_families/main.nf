#!/usr/bin/env nextflow

include { CREATE_CLUSTERS_PKL } from "${launchDir}/modules/family/main.nf"
include { REFINE_FAMILIES     } from "${launchDir}/modules/family/main.nf"

workflow GENERATE_FAMILIES {
    take:
    clusters_tsv
    refined_families_tsv
    mgnifams_fasta
    discarded_clusters

    main:
    clusters_pkl = CREATE_CLUSTERS_PKL(clusters_tsv).pkl
    refined_families = REFINE_FAMILIES(clusters_pkl, refined_families_tsv, mgnifams_fasta, discarded_clusters, params.minimum_members, params.iteration)

    emit:
    tsv = refined_families.tsv
    seed_msa = refined_families.seed_msa
    msa = refined_families.msa
    hmm = refined_families.hmm
    domtblout = refined_families.domtblout
}
