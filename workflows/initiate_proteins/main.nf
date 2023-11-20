#!/usr/bin/env nextflow

include { INITIATE_PROTEINS } from "$launchDir/subworkflows/initiate_proteins/main.nf"

workflow {
    INITIATE_PROTEINS()
}
