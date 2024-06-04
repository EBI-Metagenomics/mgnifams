#!/usr/bin/env nextflow

include { REMOVE_REDUNDANCY } from "${projectDir}/../../subworkflows/remove_redundancy/main.nf"

workflow {
    REMOVE_REDUNDANCY( params.families_output )
}
