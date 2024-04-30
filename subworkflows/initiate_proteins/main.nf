#!/usr/bin/env nextflow

include { FILTER_UNANNOTATED_SLICES } from "${params.moduleDir}/initiate.nf"
include { PUBLISH_INPUT_FASTA       } from "${params.moduleDir}/initiate.nf"

workflow INITIATE_PROTEINS {
    take:
    sequence_explorer_protein_no_header_ch

    main:
    sequence_explorer_protein_no_header_ch
        .splitText(file:true, by: params.input_csv_chunk_size)
        .set { mgy90_chunks_ch }
    fasta = FILTER_UNANNOTATED_SLICES(mgy90_chunks_ch, params.min_slice_length)
    out_fasta = fasta.collectFile(name: "mgnifams_input.fa")
    PUBLISH_INPUT_FASTA(out_fasta)

    emit:
    out_fasta
}
