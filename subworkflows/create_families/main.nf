#!/usr/bin/env nextflow

include { EXPORT_REPS; CREATE_FAMILY_FA; EXPORT_REPS_FA } from "$baseDir/modules/family.nf"
include { EXPORT_FAMILIES_CSV } from "$baseDir/modules/exporting.nf"

workflow create_families {
    take: 
    clust_tsv
    fastaFile

    main:
    reps_file = EXPORT_REPS(clust_tsv)
    EXPORT_FAMILIES_CSV(reps_file)
    reps_ch = reps_file.splitText().map { it.trim() } // removing new line chars at end of mgyps
    families = CREATE_FAMILY_FA(clust_tsv.first(), fastaFile.first(), reps_ch.take(params.debug_top_n_reps)) // TODO, remove take, take all
    known_ch = families.known_ch
    unknown_ch = families.unknown_ch
    unknown_reps_fasta = EXPORT_REPS_FA(unknown_ch.collect())

    emit:
    known_ch
    unknown_ch
    unknown_reps_fasta
}