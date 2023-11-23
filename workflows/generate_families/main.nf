#!/usr/bin/env nextflow

include { GENERATE_FAMILIES } from "$launchDir/subworkflows/generate_families/main.nf"

workflow {
    GENERATE_FAMILIES( [[id:"mmseqs_families"], params.families_tsv], params.mgnifams_input_path)
}
