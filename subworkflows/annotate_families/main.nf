#!/usr/bin/env nextflow

include { HMMSCAN } from "$baseDir/modules/hmm/hmm.nf"

workflow annotate_families {
    take:
    mafft_ch
    build_ch
    uniprot_sprot_fasta_path
    
    main:
    tblout_ch = HMMSCAN(build_ch, uniprot_sprot_fasta_path.first()).tblout_ch

    emit:
    tblout_ch
}