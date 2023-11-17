#!/usr/bin/env nextflow

include { SLICE } from "$baseDir/modules/slicing.nf"

workflow SLICE_UNANNOTATED {
    main:
    def fasta_path = params.dataDir + params.fasta_file
    def annotations_path = params.dataDir + params.annotations_file

    Channel
        .fromPath(fasta_path) 
        .set { fasta_ch }
    Channel
        .fromPath(annotations_path) 
        .set { annotations_ch }

    sliced_fasta = SLICE(fasta_ch, annotations_ch, params.min_slice_length)

    emit:
    sliced_fasta
}