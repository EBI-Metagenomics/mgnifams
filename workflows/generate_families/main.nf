#!/usr/bin/env nextflow

include { GENERATE_FAMILIES } from "${projectDir}/../../subworkflows/generate_families/main.nf"

workflow {
    GENERATE_FAMILIES( params.clusters_tsv, params.refined_families_tsv, params.mgnifams_input_fasta_path, params.discarded_clusters_txt, params.converged_families_txt, params.family_sizes_txt )
}
