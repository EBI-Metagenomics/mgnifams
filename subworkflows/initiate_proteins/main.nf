#!/usr/bin/env nextflow

include { FILTER_UNANNOTATED_SLICES } from "$launchDir/modules/initiate.nf"
include { PUBLISH_INPUT_FASTA       } from "$launchDir/modules/initiate.nf"
include { EXPORT_PROTEINS_CSV       } from "$launchDir/modules/export.nf"

workflow INITIATE_PROTEINS {
    main:
    Channel
        .fromPath(params.preprocessed_mgy90_path)
        .set { mgy90_file }

    mgy90_file
        .splitText(file:true, by: params.input_csv_chunk_size)
        .set { mgy90_chunks_ch }
    fasta = FILTER_UNANNOTATED_SLICES(mgy90_chunks_ch, params.min_slice_length)
    out_fasta = fasta.collectFile(name: "mgnifams_input.fa")
    PUBLISH_INPUT_FASTA(out_fasta)
    EXPORT_PROTEINS_CSV(out_fasta)

    emit:
    out_fasta
}
