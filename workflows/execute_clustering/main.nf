#!/usr/bin/env nextflow

include { EXECUTE_CLUSTERING } from "$launchDir/subworkflows/execute_clustering/main.nf"

workflow {
    Channel
        .fromPath(params.mgnifams_input_path) 
        .set { mgnifams_input }

    EXECUTE_CLUSTERING(mgnifams_input).families_tsv
}
