#!/usr/bin/env nextflow

// include { preprocess? } from "$baseDir/subworkflows/xxx/main.nf"
include { slice_unannotated } from "$baseDir/subworkflows/slice_unannotated/main.nf"
include { annotate_fasta } from "$baseDir/subworkflows/annotate_fasta/main.nf"

workflow {
    // preprocess?
    Channel
        .fromPath(params.msa_path) 
        .set { msa_ch }

    Channel
        .fromPath(params.tblOut_path) 
        .set { tblOut_ch }

    concat_fasta = slice_unannotated(msa_ch, tblOut_ch)
    annotate_fasta(concat_fasta)
}