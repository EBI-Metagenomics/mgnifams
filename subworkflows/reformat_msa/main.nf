#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT } from "${launchDir}/modules/hhsuite/reformat/main.nf"

workflow REFORMAT_MSA {
    take:
    fasta_ch
    
    main:
    fa_ch = HHSUITE_REFORMAT(fasta_ch, "sto", "fas").fa

    emit:
    fa_ch
}
