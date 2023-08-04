#!/usr/bin/env nextflow

include { slice_unannotated } from "$baseDir/subworkflows/slice_unannotated/main.nf"
include { annotate_slices } from "$baseDir/subworkflows/annotate_slices/main.nf"

workflow {
    sliced_fasta = slice_unannotated()
    annotate_slices(sliced_fasta)
}