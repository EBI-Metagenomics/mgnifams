#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV } from "${params.moduleDir}/export.nf"
include { QUERY_MGNPROTEIN_DB } from "${params.moduleDir}/postprocess.nf"
include { PARSE_BIOMES  } from "${params.moduleDir}/postprocess.nf"
include { PARSE_DOMAINS } from "${params.moduleDir}/postprocess.nf"

workflow {
    Channel
        .fromPath(params.updated_refined_families_path)
        .set { updated_refined_families_ch }

    EXPORT_MGNIFAMS_CSV(params.mgnifams_output_dir)

    query_results = QUERY_MGNPROTEIN_DB(params.db_config_file, updated_refined_families_ch)

    PARSE_BIOMES(query_results)
    PARSE_DOMAINS(query_results, params.updated_refined_families_path)
}
