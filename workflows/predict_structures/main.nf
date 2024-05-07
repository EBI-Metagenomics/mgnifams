#!/usr/bin/env nextflow

include { PREDICT_STRUCTURES } from "${projectDir}/../../subworkflows/predict_structures/main.nf"

workflow {
    Channel
        .fromPath(params.msa_sto_path) 
        .set { msa_sto_ch }

    PREDICT_STRUCTURES(msa_sto_ch)
}
