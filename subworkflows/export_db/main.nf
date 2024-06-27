#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV } from "${params.moduleDir}/export.nf"
include { QUERY_MGNPROTEIN_DB } from "${params.moduleDir}/postprocess.nf"
include { PARSE_BIOMES        } from "${params.moduleDir}/postprocess.nf"
include { PARSE_DOMAINS       } from "${params.moduleDir}/postprocess.nf"
// TODO
include { INITIATE_SQLITE     } from "${params.moduleDir}/postprocess.nf"
// include { UPDATE_BLOBS_PARALLEL } from "${params.moduleDir}/postprocess.nf"

workflow EXPORT_DB {
    take:
    fam_metadata
    fam_converged
    refined_families
    pfam_hits
    foldseek_hits
    predict_scores
    
    main:
    tables = EXPORT_MGNIFAMS_CSV(fam_metadata, fam_converged, \
        refined_families, pfam_hits, foldseek_hits, predict_scores)
    query_results  = QUERY_MGNPROTEIN_DB(params.db_config_file, refined_families)
    biome_results  = PARSE_BIOMES(query_results)
    domain_results = PARSE_DOMAINS(query_results, refined_families)

    INITIATE_SQLITE(params.db_schema_file, tables)
}
