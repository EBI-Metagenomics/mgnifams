#!/usr/bin/env nextflow

include { INITIATE_PROTEINS } from "$launchDir/subworkflows/initiate_proteins/main.nf"

workflow {
    Channel
        .fromPath(params.preprocessed_mgy90_path)
        .set { preprocessed_mgy90 }

    INITIATE_PROTEINS(preprocessed_mgy90)
}
