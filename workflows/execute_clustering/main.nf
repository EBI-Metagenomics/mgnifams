#!/usr/bin/env nextflow

include { EXECUTE_CLUSTERING } from "$launchDir/subworkflows/execute_clustering/main.nf"

workflow {
    EXECUTE_CLUSTERING(params.mgnifams_input_path).families_tsv
}
