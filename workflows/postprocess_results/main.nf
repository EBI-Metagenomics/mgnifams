#!/usr/bin/env nextflow

include { PARSE_BIOMES  } from "${params.moduleDir}/postprocess.nf"
include { PARSE_DOMAINS } from "${params.moduleDir}/postprocess.nf"

workflow {
    PARSE_BIOMES(params.query_results_path)
    PARSE_DOMAINS(params.query_results_path, params.refined_families_path)
}
