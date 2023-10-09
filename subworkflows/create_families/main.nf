#!/usr/bin/env nextflow

include { PARSE_FAMILIES } from "$baseDir/modules/family.nf"
// include { EXPORT_FAMILIES_CSV } from "$baseDir/modules/exporting.nf"

workflow create_families {
    take:
    fastaFile
    clust_tsv
    
    main:
    families = PARSE_FAMILIES(fastaFile, clust_tsv)
    reps_ids = families.reps_ids
    reps_fasta = families.reps_fasta
    singleton_ids = families.singleton_ids
    non_singletons_ch = families.non_singletons_ch
    // EXPORT_FAMILIES_CSV(unknown_ids_file, 'unknown')

    emit:
    reps_ids
    reps_fasta
    singleton_ids
    non_singletons_ch
}