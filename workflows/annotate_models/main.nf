#!/usr/bin/env nextflow

include { ANNOTATE_MODELS } from "$launchDir/subworkflows/annotate_models/main.nf"

workflow {
    Channel
        .fromPath(params.msa_folder_path)
        .map { filepath ->
            return [ [id:"annotated_models"], file(filepath) ]
        }
        .set { msa_ch }
    
    ANNOTATE_MODELS( msa_ch, params.hhdb_folder_path )
}
