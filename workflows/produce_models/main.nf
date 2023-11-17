#!/usr/bin/env nextflow

include { PRODUCE_MODELS } from "../../subworkflows/produce_models/main.nf"

workflow {
    Channel
        .fromPath(params.families_path) 
        .set { families_ch }
    
    PRODUCE_MODELS(families_ch)
}
