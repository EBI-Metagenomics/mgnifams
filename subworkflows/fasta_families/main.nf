#!/usr/bin/env nextflow

include { EXPORT_REPS, CREATE_FAMILY_FA } from "$baseDir/modules/families/family.nf"

workflow families {
    take: 
        clust_tsv,
        fastaFile

    main:
        CREATE_FAMILY_FA(clust_tsv, fastaFile) // todo export reps first, then redo

    emit:
        CREATE_FAMILY_FA.out
}