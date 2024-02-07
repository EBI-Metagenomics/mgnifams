#!/usr/bin/env nextflow

include { REFORMAT_MSA as REFORMAT_SEED_MSA    } from "$launchDir/subworkflows/reformat_msa/main.nf"
include { REFORMAT_MSA as REFORMAT_HMMALIGN_MSA} from "$launchDir/subworkflows/reformat_msa/main.nf"

workflow {
    Channel
        .fromPath(params.seed_msa_folder_path)
        .map { filepath ->
            return [ [id:"seed_msa"], file(filepath) ]
        }
        .set { seed_msa_ch }
    
    Channel
        .fromPath(params.msa_folder_path)
        .map { filepath ->
            return [ [id:"msa"], file(filepath) ]
        }
        .set { hmmalign_msa_ch }
    
    REFORMAT_SEED_MSA(seed_msa_ch)
    REFORMAT_HMMALIGN_MSA( hmmalign_msa_ch )
}
