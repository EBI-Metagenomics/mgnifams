#!/usr/bin/env nextflow

include { SLICE } from "$baseDir/modules/slicing.nf"

workflow slice_unannotated {
    take:
    fasta_path
    annotations_path
    min_slice_length
    
    main:
    sliced_fasta = SLICE(fasta_path, annotations_path, min_slice_length)

    emit:
    sliced_fasta
}