#!/usr/bin/env nextflow

include { PREPROCESS_INPUT   } from "$launchDir/subworkflows/preprocess_input/main.nf"
include { INITIATE_PROTEINS  } from "$launchDir/subworkflows/initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING } from "$launchDir/subworkflows/execute_clustering/main.nf"

workflow {
    Channel
        .fromPath(params.mgy90_path)
        .set { mgy90_file_bz2 }

    preprocessed_mgy90_file = PREPROCESS_INPUT(mgy90_file_bz2).preprocessed_mgy90_file
    out_fasta               = INITIATE_PROTEINS(preprocessed_mgy90_file).out_fasta
    families_tsv            = EXECUTE_CLUSTERING(out_fasta).families_tsv
}
