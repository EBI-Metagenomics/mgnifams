#!/usr/bin/env nextflow

include { EXPORT_REPS; CREATE_FAMILY_FA } from "$baseDir/modules/families/family.nf"

workflow fasta_families {
    take: 
    clust_tsv
    fastaFile

    main: // TODO
    reps_ch = EXPORT_REPS(clust_tsv)
    // CREATE_FAMILY_FA(reps_ch, fastaFile)

    // emit:
    // CREATE_FAMILY_FA.out
}