#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV } from "${params.moduleDir}/export.nf"
include { QUERY_MGNPROTEIN_DB } from "${params.moduleDir}/postprocess.nf"
include { PARSE_BIOMES        } from "${params.moduleDir}/postprocess.nf"
include { PARSE_DOMAINS       } from "${params.moduleDir}/postprocess.nf"

workflow EXPORT_DATA {
    take:
    fam_metadata
    fam_converged
    refined_families
    pfam_hits
    foldseek_hits
    predict_scores
    
    main:
    EXPORT_MGNIFAMS_CSV( fam_metadata, fam_converged, \
        refined_families, pfam_hits, foldseek_hits, predict_scores )
    query_results = QUERY_MGNPROTEIN_DB(params.db_config_file, refined_families)
    PARSE_BIOMES(query_results)
    PARSE_DOMAINS(query_results, refined_families)
}
