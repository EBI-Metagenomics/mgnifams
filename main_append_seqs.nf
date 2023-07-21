#!/usr/bin/env nextflow

// include { execute_clustering } from "$baseDir/subworkflows/execute_clustering/main.nf"
// include { create_families } from "$baseDir/subworkflows/create_families/main.nf"
// include { produce_models } from "$baseDir/subworkflows/produce_models/main.nf"
// include { annotate_families } from "$baseDir/subworkflows/annotate_families/main.nf"

workflow {
    // Channel
    //     .fromPath(params.fasta_path) 
    //     .set { fastaFile }

    // mmseqs = execute_clustering(fastaFile)
    // families_ch = create_families(mmseqs.clu_tsv, fastaFile)
    // models = produce_models(families_ch)
    // annotate_families(mmseqs.rep_fa, models.build_ch)
}