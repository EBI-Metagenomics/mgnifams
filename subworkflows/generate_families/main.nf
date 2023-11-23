#!/usr/bin/env nextflow

include { REFINE_FAMILIES } from "${launchDir}/modules/family.nf"

workflow GENERATE_FAMILIES {
    take:
    families_tsv
    mgnifams_fasta
    
    main:
    refined_families = REFINE_FAMILIES(families_tsv, mgnifams_fasta, params.minimum_members).tsv

    emit:
    refined_families
}
