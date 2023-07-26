#!/usr/bin/env nextflow

include { CONCAT_FILES } from "$baseDir/modules/general.nf"
include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
include { produce_models } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_families } from "$baseDir/subworkflows/annotate_families/main.nf"

workflow {
    Channel
        .fromPath(params.mgnify_fasta_path) 
        .set { mgnify_fasta_file }
    Channel
        .fromPath(params.uniprot_sprot_fasta_input_debug_path) 
        .set { uniprot_sprot_fasta_input_debug_file } // TODO change to main Uniprot file

    combined_fasta_file = CONCAT_FILES(mgnify_fasta_file, uniprot_sprot_fasta_input_debug_file)
    mmseqs = execute_clustering(combined_fasta_file)
    families_ch = create_families(mmseqs.clu_tsv, combined_fasta_file)
    models = produce_models(families_ch)
    annotate_families(mmseqs.rep_fa, models.build_ch)
}