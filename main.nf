#!/usr/bin/env nextflow

include { initiate_proteins } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
include { produce_models as produce_unknown_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { produce_models as produce_known_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_unknown } from "$baseDir/subworkflows/annotate_unknown/main.nf"
include { annotate_known } from "$baseDir/subworkflows/annotate_known/main.nf"

workflow {
    combined_fasta_file = initiate_proteins()
    mmseqs = execute_clustering(combined_fasta_file)
    families = create_families(mmseqs.clu_tsv, combined_fasta_file)
    // unknown_models = produce_unknown_models(families.unknown_ch)
    // annotate_unknown(families.unknown_reps_fasta, unknown_models.build_ch)
}