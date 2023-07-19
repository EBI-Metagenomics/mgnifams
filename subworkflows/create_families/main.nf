#!/usr/bin/env nextflow

include { EXPORT_REPS; CREATE_FAMILY_FA;
          EXPORT_REPS_FA } from "$baseDir/modules/families/family.nf"

workflow create_families {
    take: 
    clust_tsv
    fastaFile

    main:
    reps_ch = EXPORT_REPS(clust_tsv).splitText().map { it.trim() } // removing new line chars at end of mgyps
    families_ch = CREATE_FAMILY_FA(clust_tsv.first(), fastaFile.first(), reps_ch.take(20)) // TODO, remove top 20
    reps_fasta = EXPORT_REPS_FA(families_ch.collect())

    emit:
    families_ch
    reps_fasta
}