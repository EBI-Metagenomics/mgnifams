#!/usr/bin/env nextflow

include { PREDICT_STRUCTURES  } from "${projectDir}/../../subworkflows/predict_structures/main.nf"
include { ANNOTATE_STRUCTURES } from "${projectDir}/../../subworkflows/annotate_structures/main.nf"

workflow {
    Channel
        .fromPath(params.msa_sto_path) 
        .map { filepath ->
            return [ [id:"msa_sto"], file(filepath) ]
        }
        .set { msa_sto_ch }

    pdb_ch = PREDICT_STRUCTURES(msa_sto_ch).pdb_ch
    ANNOTATE_STRUCTURES(pdb_ch)
}
