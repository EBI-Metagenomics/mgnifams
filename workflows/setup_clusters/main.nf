#!/usr/bin/env nextflow

include { PREPROCESS_INPUT   } from "${projectDir}/../../subworkflows/preprocess_input/main.nf"
include { INITIATE_PROTEINS  } from "${projectDir}/../../subworkflows/initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING } from "${projectDir}/../../subworkflows/execute_clustering/main.nf"

workflow {
    preprocessed_sequence_explorer_protein_ch = PREPROCESS_INPUT(params.sequence_explorer_protein_path, params.compress_mode).preprocessed_sequence_explorer_protein_ch
    fasta_ch = INITIATE_PROTEINS( preprocessed_sequence_explorer_protein_ch ).fasta_ch
    EXECUTE_CLUSTERING( fasta_ch )
}
