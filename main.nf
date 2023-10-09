#!/usr/bin/env nextflow

include { initiate_proteins } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
include { produce_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_slices } from "$baseDir/subworkflows/annotate_slices/main.nf"
include { FIND_UNANNOTATED_IDS } from "$baseDir/modules/general.nf"

workflow {
    combined_fasta_file = initiate_proteins()
    mmseqs = execute_clustering(combined_fasta_file)
    families = create_families(combined_fasta_file, mmseqs.clu_tsv)
    unknown_models = produce_models(families.non_singletons_ch)
    // annotations_ch = annotate_slices(families.unknown_reps_fasta)
    // FIND_UNANNOTATED_IDS(annotations_ch, families.reps_file)
}