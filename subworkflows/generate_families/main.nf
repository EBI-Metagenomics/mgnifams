#!/usr/bin/env nextflow

include { REFINE_FAMILIES } from "${launchDir}/modules/family/main.nf"

workflow GENERATE_FAMILIES {
    take:
    families_tsv
    mgnifams_fasta
    // families_pkl
    // mgnifams_pkl
    
    main:
    refined_families = REFINE_FAMILIES(families_tsv, mgnifams_fasta, params.minimum_members)
    // refined_families = REFINE_FAMILIES(families_tsv, mgnifams_fasta, params.minimum_members, families_pkl, mgnifams_pkl) // from saved state

    emit:
    tsv = refined_families.tsv
    seed_msa = refined_families.seed_msa
    msa = refined_families.msa
    hmm = refined_families.hmm
    domtblout = refined_families.domtblout
    // log = refined_families.log
}
