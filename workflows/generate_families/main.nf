#!/usr/bin/env nextflow

include { EXECUTE_CLUSTERING } from "$launchDir/subworkflows/execute_clustering/main.nf"
// include { CREATE_FAMILIES    } from "$launchDir/subworkflows/create_families/main.nf"

workflow {
    families_tsv = EXECUTE_CLUSTERING(params.mgnifams_input_path).families_tsv
    // CREATE_FAMILIES(families_tsv)
}
