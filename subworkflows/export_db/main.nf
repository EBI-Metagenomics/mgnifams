#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV        } from "${params.moduleDir}/export.nf"
include { INITIATE_SQLITE            } from "${params.moduleDir}/postprocess.nf"
include { QUERY_MGNPROTEIN_DB        } from "${params.moduleDir}/postprocess.nf"
include { PARSE_BIOMES               } from "${params.moduleDir}/postprocess.nf"
include { CHUNK_QUERY_RESULTS_FOLDER } from "${params.moduleDir}/postprocess.nf"
include { PARSE_DOMAINS              } from "${params.moduleDir}/postprocess.nf"
include { APPEND_BLOBS               } from "${params.moduleDir}/postprocess.nf"

workflow EXPORT_DB {
    take:
    fam_metadata
    fam_converged
    refined_families
    pfam_hits
    foldseek_hits
    predict_scores
    cif_ch
    seed_msa_ch
    msa_ch
    hmm_ch
    rf_ch
    
    main:
    tables = EXPORT_MGNIFAMS_CSV(fam_metadata, fam_converged, \
        refined_families, pfam_hits, foldseek_hits, predict_scores)
    db = INITIATE_SQLITE(params.db_schema_file, tables)

    query_results  = QUERY_MGNPROTEIN_DB(params.db_config_file, refined_families)
    biome_results  = PARSE_BIOMES(query_results)
    chunks         = CHUNK_QUERY_RESULTS_FOLDER(query_results, params.query_results_chunk_size)
    domain_results = PARSE_DOMAINS(chunks.flatten(), query_results.first(), refined_families.first())
    APPEND_BLOBS(db.first(), cif_ch, seed_msa_ch.first(), msa_ch.first(), hmm_ch.first(), rf_ch.first(), \
        biome_results.first(), domain_results.flatten(), params.update_blobs_from_row_id)
}
