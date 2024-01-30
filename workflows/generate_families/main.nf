#!/usr/bin/env nextflow

include { GENERATE_FAMILIES } from "$launchDir/subworkflows/generate_families/main.nf"

workflow {
    GENERATE_FAMILIES( params.clusters_tsv, params.refined_families_tsv, params.mgnifams_dict_fasta_path, params.discarded_clusters_txt )
}
