#!/usr/bin/env nextflow

include { CREATE_DB; LINCLUST; CREATE_TSV } from "$baseDir/modules/mmseqs2.nf"

workflow execute_clustering {
    take:
    fastaFile

    main:
    db_ch = CREATE_DB(fastaFile)
    clust_ch = LINCLUST(db_ch)
    CREATE_TSV(db_ch, clust_ch)
        
    emit:
    CREATE_TSV.out
}