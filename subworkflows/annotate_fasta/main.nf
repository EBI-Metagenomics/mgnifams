#!/usr/bin/env nextflow

include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
// TODO more annotation sources: Uniprot, eggNOG, protENN

workflow annotate_fasta {
    take:
    fasta_ch
    
    main:
    INTERPROSCAN(fasta_ch)

    // emit:
}