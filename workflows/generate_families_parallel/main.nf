#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL } from "${projectDir}/../../subworkflows/generate_families_parallel/main.nf"

workflow {
    GENERATE_FAMILIES_PARALLEL( params.clusters_tsv, params.checked_clusters_txt, params.mgnifams_input_fasta_path)
}
