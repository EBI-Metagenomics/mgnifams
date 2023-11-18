#!/usr/bin/env nextflow

include { FILTER_FAMILIES  } from "$baseDir/modules/family.nf"
include { PARSE_FAMILIES } from "$baseDir/modules/family.nf"

workflow CREATE_FAMILIES {
    take:
    fastaFile
    clust_tsv
    
    main:
    filtered_families = FILTER_FAMILIES(clust_tsv, params.family_member_threshold)
    families = PARSE_FAMILIES(fastaFile, filtered_families.filtered_clusters)
    reps_ids = families.reps_ids
    reps_fasta = families.reps_fasta
    families_folder = families.families_folder

    emit:
    reps_ids
    reps_fasta
    families_folder
}
