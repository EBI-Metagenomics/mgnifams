#!/usr/bin/env nextflow

include { INITIATE_PROTEINS } from "${projectDir}/../../subworkflows/initiate_proteins/main.nf"

workflow {
    Channel
        .fromPath(params.sequence_explorer_protein_no_header_path)
        .set { sequence_explorer_protein_no_header_ch }

    INITIATE_PROTEINS(sequence_explorer_protein_no_header_ch)
}
