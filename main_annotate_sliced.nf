#!/usr/bin/env nextflow

// include { preprocess? } from "$baseDir/subworkflows/xxx/main.nf"
include { slice_unannotated } from "$baseDir/subworkflows/slice_unannotated/main.nf"
// include { annotate_families } from "$baseDir/subworkflows/annotate_families/main.nf"

workflow {
    // preprocess?
    Channel
        .fromPath(params.msa_path) 
        .set { msa_ch }

    Channel
        .fromPath(params.tblOut_path) 
        .set { tblOut_ch }

    slice_unannotated(msa_ch, tblOut_ch)
    // annotate_families(models.mafft_ch, models.build_ch, uniprot_sprot_fasta_path)
}