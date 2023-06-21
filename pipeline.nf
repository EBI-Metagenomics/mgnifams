#!/usr/bin/env nextflow

params.fasta_path = "data/mgnify_500K_proteins.fa.gz" // "data/cluster_reps_50k.fa"

include { CREATE_DB; LINCLUST; CREATE_TSV } from './modules/mmseqs2/mmseqs2.nf'

workflow {
    Channel
        .fromPath(params.fasta_path) 
        .set { fastaFile }

    db_ch = CREATE_DB(fastaFile)
    clust_ch = LINCLUST(db_ch)
    CREATE_TSV(db_ch, clust_ch)
}