#!/usr/bin/env nextflow

include { QUERY_MGNPROTEIN_DB } from "${params.moduleDir}/postprocess.nf"

workflow {
    Channel
        .fromPath(params.updated_refined_families_path)
        .set { updated_refined_families_ch }

    QUERY_MGNPROTEIN_DB(params.db_config_file, updated_refined_families_ch)
}
