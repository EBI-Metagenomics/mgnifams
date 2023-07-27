#!/usr/bin/env nextflow

include { initiate_proteins } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
include { produce_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_families } from "$baseDir/subworkflows/annotate_families/main.nf"

workflow {
    combined_fasta_file = initiate_proteins()
    mmseqs = execute_clustering(combined_fasta_file)
    families_all = create_families(mmseqs.clu_tsv, combined_fasta_file)
    models = produce_models(families_all)
    // TODO separate known and unknown fmailies after this point
    annotate_families(mmseqs.rep_fa, models.build_ch)
}