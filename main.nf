#!/usr/bin/env nextflow

include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
include { produce_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_families } from "$baseDir/subworkflows/annotate_families/main.nf"

workflow {
    Channel
        .fromPath(params.fasta_path) 
        .set { fastaFile }

    Channel
        .fromPath(params.uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_path }

    cluster_tsv_ch = execute_clustering(fastaFile)
    families_ch = create_families(cluster_tsv_ch, fastaFile)
    models = produce_models(families_ch)
    annotate_families(models.mafft_ch, models.build_ch, uniprot_sprot_fasta_path)
}