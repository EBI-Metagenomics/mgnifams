#!/usr/bin/env nextflow

include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"

workflow annotate_slices {
    take:
    fasta_ch
    
    main:
    interpro_ch = INTERPROSCAN(fasta_ch)
}