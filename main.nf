#!/usr/bin/env nextflow

include { initiate_proteins } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
include { produce_models as produce_unknown_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { produce_models as produce_known_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_families } from "$baseDir/subworkflows/annotate_families/main.nf"

workflow {
    combined_fasta_file = initiate_proteins()
    mmseqs = execute_clustering(combined_fasta_file)
    families = create_families(mmseqs.clu_tsv, combined_fasta_file)
    unknown_models = produce_unknown_models(families.unknown_ch)
    known_models = produce_known_models(families.known_ch)
    // TODO separate known and unknown families annotation
    annotate_families(families.unknown_reps_fasta, unknown_models.build_ch)
}