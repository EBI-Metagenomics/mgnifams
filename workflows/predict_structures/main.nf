#!/usr/bin/env nextflow

include { PREDICT_STRUCTURES } from "${launchDir}/subworkflows/predict_structures/main.nf"

workflow {
    Channel
        .fromPath(params.msa_path) 
        .set { msa }
        
    Channel
        .fromPath(params.wp1_unannotated_ids_path) 
        .set { wp1_unannotated_ids }

    PREDICT_STRUCTURES(msa, wp1_unannotated_ids, "all")
}
