#!/usr/bin/env nextflow

include { REMOVE_REDUNDANCY } from "${projectDir}/../../subworkflows/remove_redundancy/from_folder.nf"

workflow {
    Channel
        .fromPath(params.seed_msa_sto_path)
        .map { filepath ->
            return [ [id:"remove_redundancy"], file(filepath) ]
        }
        .set { seed_msa_sto_ch }

    REMOVE_REDUNDANCY( seed_msa_sto_ch, params.families_output )
}
