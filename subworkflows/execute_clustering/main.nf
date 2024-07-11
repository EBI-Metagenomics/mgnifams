#!/usr/bin/env nextflow

include { MMSEQS_CREATEDB     } from "../../modules/mmseqs/createdb/main.nf"
include { MMSEQS_LINCLUST     } from "../../modules/mmseqs/linclust/main.nf"
include { EXPORT_CLUSTERS_TSV } from "../../modules/mmseqs/export_clusters_tsv.nf"

workflow EXECUTE_CLUSTERING {
    take:
    mgnifams_input_fasta

    main:
    mmseqs_db       = MMSEQS_CREATEDB( mgnifams_input_fasta.map { fasta -> [ [ id:'mmseqs_db' ], fasta ] } ).db
    mmseqs_clusters = MMSEQS_LINCLUST(mmseqs_db).db_cluster
    clusters        = EXPORT_CLUSTERS_TSV(mmseqs_db, mmseqs_clusters)
    
    emit:
    clusters_tsv  = clusters.tsv
    num_sequences = clusters.num_sequences
}
