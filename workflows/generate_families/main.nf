#!/usr/bin/env nextflow

include { GENERATE_FAMILIES } from "$launchDir/subworkflows/generate_families/main.nf"

workflow {
    GENERATE_FAMILIES( params.families_tsv, params.mgnifams_dict_fasta_path, params.mode )
}
