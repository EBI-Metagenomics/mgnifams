#!/usr/bin/env nextflow

include { EXPORT_REPS; CREATE_FAMILY_FA; KEEP_UNKNOWN; KEEP_KNOWN } from "$baseDir/modules/family.nf"
include { EXPORT_FAMILIES_CSV } from "$baseDir/modules/exporting.nf"

workflow create_families {
    take: 
    clust_tsv
    fastaFile

    main:
    reps_file = EXPORT_REPS(clust_tsv)
    EXPORT_FAMILIES_CSV(reps_file)
    reps_ch = reps_file.splitText().map { it.trim() } // removing new line chars at end of mgyps
    families_all = CREATE_FAMILY_FA(clust_tsv.first(), fastaFile.first(), reps_ch.take(params.debug_top_n_reps)) // TODO, remove take, take all
    // families_ch = KEEP_UNKNOWN(families_all)
    // KEEP_KNOWN(families_all)

    emit:
    families_all
}