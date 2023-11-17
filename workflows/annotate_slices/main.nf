#!/usr/bin/env nextflow

include { SLICE_UNANNOTATED } from "$baseDir/subworkflows/slice_unannotated/main.nf"
include { ANNOTATE_SLICES } from "$baseDir/subworkflows/annotate_slices/main.nf"

workflow {
    sliced_fasta = SLICE_UNANNOTATED()
    ANNOTATE_SLICES(sliced_fasta)
}