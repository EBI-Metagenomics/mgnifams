#!/usr/bin/env nextflow

include { FILTER_UNANNOTATED_SLICES } from "${params.moduleDir}/initiate.nf"
include { PUBLISH_INPUT_FASTA       } from "${params.moduleDir}/initiate.nf"

workflow INITIATE_PROTEINS {
    take:
    sequence_explorer_protein_no_header_ch

    main:
    sequence_explorer_protein_no_header_ch
        .splitText(file:true, by: params.input_csv_chunk_size)
        .set { sequence_chunk_ch }
    fasta_chunk_ch = FILTER_UNANNOTATED_SLICES(sequence_chunk_ch, params.min_sequence_length)
    fasta_ch = fasta_chunk_ch.collectFile(name: "mgnifams_input.fa")
    PUBLISH_INPUT_FASTA(fasta_ch)

    emit:
    fasta_ch
}
