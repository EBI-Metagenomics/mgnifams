#!/usr/bin/env nextflow

include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"

workflow annotate_slices {
    take:
    fasta_ch
    
    main:
    Channel
        .fromPath(params.eggnong_data_dir) 
        .set { ch_eggnong_data_dir }
    Channel
        .fromPath(params.eggnong_diamond_db) 
        .set { ch_eggnog_diamond_db }
    Channel
        .fromPath(params.eggnog_db) 
        .set { ch_eggnog_db }

    interpro_ch = INTERPROSCAN(fasta_ch)
    eggnog_ch = EGGNOG_MAPPER(fasta_ch, ch_eggnong_data_dir, ch_eggnog_diamond_db, ch_eggnog_db)
}