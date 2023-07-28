#!/usr/bin/env nextflow

include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
// TODO more annotation sources: Uniprot, eggNOG, protENN, probably merge part of it with annotate_family subworkflow (instead of hmmscan Uniprot)

workflow annotate_fasta {
    take:
    fasta_ch
    
    main:
    INTERPROSCAN(fasta_ch)

    // emit:
}