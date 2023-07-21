#!/usr/bin/env nextflow

include { update_clusters } from "$baseDir/subworkflows/update_clusters/main.nf"

workflow {
    Channel
        .fromPath(params.new_fasta_path) 
        .set { new_fasta_file }

    Channel
        .fromPath(params.db_path) 
        .set { old_db_ch }

    Channel
        .fromPath(params.clu_path) 
        .set { old_clu_ch }

    update_clusters(new_fasta_file, old_db_ch.collect(), old_clu_ch.collect())
}