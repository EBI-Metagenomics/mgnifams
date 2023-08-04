#!/usr/bin/env nextflow

include { slice_unannotated } from "$baseDir/subworkflows/slice_unannotated/main.nf"
include { annotate_slices } from "$baseDir/subworkflows/annotate_slices/main.nf"

workflow {
    Channel
        .fromPath(params.fasta_file) 
        .set { fasta_path }
    Channel
        .fromPath(params.annotations_file) 
        .set { annotations_path }

    sliced_fasta = slice_unannotated(fasta_path, annotations_path, params.min_slice_length)
    // annotate_slices(sliced_fasta)
}