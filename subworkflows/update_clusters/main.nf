#!/usr/bin/env nextflow

include { CREATEDB; CLUSTERUPDATE } from "$baseDir/modules/mmseqs2.nf"

workflow UPDATE_CLUSTERS {
    take:
    new_fasta_ch
    old_db_ch
    old_clu_ch
    
    main:
    new_db_ch = CREATEDB(new_fasta_ch)
    updated_clu_ch = CLUSTERUPDATE(old_db_ch, new_db_ch, old_clu_ch)
}