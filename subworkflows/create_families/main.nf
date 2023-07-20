#!/usr/bin/env nextflow

include { EXPORT_REPS; CREATE_FAMILY_FA } from "$baseDir/modules/family.nf"

workflow create_families {
    take: 
    clust_tsv
    fastaFile

    main:
    reps_ch = EXPORT_REPS(clust_tsv).splitText().map { it.trim() } // removing new line chars at end of mgyps
    families_ch = CREATE_FAMILY_FA(clust_tsv.first(), fastaFile.first(), reps_ch.take(params.debug_top_n_reps)) // TODO, remove take, take all

    emit:
    families_ch
}