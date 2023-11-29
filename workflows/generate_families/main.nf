#!/usr/bin/env nextflow

include { GENERATE_FAMILIES } from "$launchDir/subworkflows/generate_families/main.nf"

workflow {
    // GENERATE_FAMILIES( params.families_tsv, params.mgnifams_input_path )
    GENERATE_FAMILIES( params.families_tsv, params.mgnifams_input_path, params.families_pkl )
    // GENERATE_FAMILIES( params.families_tsv, params.mgnifams_input_path, params.families_pkl, params.mgnifams_input_pkl )
}
