#!/usr/bin/env nextflow

include { REFORMAT_MSA } from "$launchDir/subworkflows/reformat_msa/main.nf"

workflow {
    Channel
        .fromPath(params.msa_folder_path)
        .map { filepath ->
            return [ [id:"reformat_msa"], file(filepath) ]
        }
        .set { msa_ch }
    
    REFORMAT_MSA( msa_ch )
}
