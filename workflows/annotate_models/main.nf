#!/usr/bin/env nextflow

include { ANNOTATE_MODELS } from "$launchDir/subworkflows/annotate_models/main.nf"

workflow {
    Channel
        .fromPath(params.seed_msa_path)
        .map { filepath ->
            return [ [id:"annotated_models"], file(filepath) ]
        }
        .set { seed_msa_ch }
    
    ANNOTATE_MODELS( seed_msa_ch, "hhblits" )
}
