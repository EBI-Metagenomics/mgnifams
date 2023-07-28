#!/usr/bin/env nextflow

include { CREATEDB; LINCLUST; CREATETSV } from "$baseDir/modules/mmseqs2.nf"
include { EXPORT_CLUSTERING_CSV } from "$baseDir/modules/exporting.nf"

workflow execute_clustering {
    take:
    fasta_file

    main:
    db_ch = CREATEDB(fasta_file)
    clu_ch = LINCLUST(db_ch)
    clu_tsv = CREATETSV(db_ch, clu_ch)
    EXPORT_CLUSTERING_CSV(clu_tsv)
    // rep_fa = CONVERT2FASTA(db_ch, clu_ch)
    
    emit:
    clu_tsv
    // rep_fa
}