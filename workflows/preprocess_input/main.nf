#!/usr/bin/env nextflow

include { PREPROCESS_INPUT } from "$launchDir/subworkflows/preprocess_input/main.nf"

workflow {
    PREPROCESS_INPUT()
}
