#!/usr/bin/env nextflow

include { CREATEDB; LINCLUST; CREATETSV; CONVERT2FASTA } from "$baseDir/modules/mmseqs2.nf"

workflow execute_clustering {
    take:
    fastaFile

    main:
    db_ch = CREATEDB(fastaFile)
    clu_ch = LINCLUST(db_ch)
    clu_tsv = CREATETSV(db_ch, clu_ch)
    rep_fa = CONVERT2FASTA(db_ch, clu_ch)
    
    emit:
    clu_tsv
    rep_fa
}