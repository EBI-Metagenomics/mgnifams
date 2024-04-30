#!/usr/bin/env nextflow

include { PREPROCESS_INPUT } from "${projectDir}/../../subworkflows/preprocess_input/main.nf"

workflow {
    Channel
        .fromPath(params.sequence_explorer_protein_path)
        .set { sequence_explorer_protein_ch }

    PREPROCESS_INPUT(sequence_explorer_protein_ch, params.compress_mode)
}
