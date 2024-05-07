#!/usr/bin/env nextflow

include { ANNOTATE_MODELS } from "${projectDir}/../../subworkflows/annotate_models/main.nf"

workflow {
    Channel
        .fromPath(params.seed_msa_path)
        .map { filepath ->
            return [ [id:"annotate_models"], file(filepath) ]
        }
        .set { seed_msa_ch }
    
    ANNOTATE_MODELS( seed_msa_ch, params.hh_mode )
}
