#!/usr/bin/env nextflow

include { KEEP_FAMILIES  } from "$baseDir/modules/family.nf"
include { PARSE_FAMILIES } from "$baseDir/modules/family.nf"
// include { EXPORT_FAMILIES_CSV } from "$baseDir/modules/exporting.nf"

workflow create_families {
    take:
    fastaFile
    clust_tsv
    
    main:
    filtered_families = KEEP_FAMILIES(clust_tsv, params.family_member_threshold)
    families = PARSE_FAMILIES(fastaFile, filtered_families.filtered_clusters)
    reps_ids = families.reps_ids
    reps_fasta = families.reps_fasta
    families_folder = families.families_folder
    // EXPORT_FAMILIES_CSV(unknown_ids_file, 'unknown')

    emit:
    reps_ids
    reps_fasta
    families_folder
}