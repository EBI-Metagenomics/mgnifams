#!/usr/bin/env nextflow

include { HMMSCAN } from "$baseDir/modules/hmm.nf"
include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"
// TODO more annotation sources: eggNOG, protENN

workflow annotate_families {
    take:
    reps_fa
    build_ch
    
    main:
    Channel
        .fromPath(params.uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_path }

    INTERPROSCAN(reps_fa)
    EGGNOG_MAPPER(reps_fa)
    tblout_ch = HMMSCAN(build_ch, uniprot_sprot_fasta_path.first()).tblout_ch

    emit:
    tblout_ch
}