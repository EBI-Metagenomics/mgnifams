#!/usr/bin/env nextflow

include { PREPROCESS_INPUT } from "$launchDir/subworkflows/preprocess_input/main.nf"

workflow {
    Channel
        .fromPath(params.mgy90_path)
        .set { mgy90_file_bz2 }

    PREPROCESS_INPUT(mgy90_file_bz2)
}
