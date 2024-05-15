#!/usr/bin/env nextflow

include { REFORMAT_MSA as REFORMAT_SEED_MSA     } from "${projectDir}/../../subworkflows/reformat_msa/main.nf"
include { REFORMAT_MSA as REFORMAT_HMMALIGN_MSA } from "${projectDir}/../../subworkflows/reformat_msa/main.nf"
include { ANNOTATE_MODELS                       } from "${projectDir}/../../subworkflows/annotate_models/main.nf"
include { PREDICT_STRUCTURES                    } from "${projectDir}/../../subworkflows/predict_structures/main.nf"
include { ANNOTATE_STRUCTURES                   } from "${projectDir}/../../subworkflows/annotate_structures/main.nf"

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
    
    fa_seed_msa_ch = REFORMAT_SEED_MSA(seed_msa_ch).fa_ch
    REFORMAT_HMMALIGN_MSA( hmmalign_msa_ch )
    ANNOTATE_MODELS( fa_seed_msa_ch, params.hh_mode )

    pdb_ch = PREDICT_STRUCTURES(hmmalign_msa_ch).pdb_ch
    ANNOTATE_STRUCTURES(pdb_ch)
}
