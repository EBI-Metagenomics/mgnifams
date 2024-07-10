#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV   } from "../../modules/export.nf"
include { QUERY_MGNPROTEIN_DB   } from "../../modules/postprocess.nf"
include { PARSE_BIOMES          } from "../../modules/postprocess.nf"
include { PARSE_DOMAINS         } from "../../modules/postprocess.nf"
include { INITIATE_SQLITE       } from "../../modules/postprocess.nf"
include { APPEND_BLOBS_PARALLEL } from "../../modules/postprocess.nf"

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
    query_results  = QUERY_MGNPROTEIN_DB(params.db_config_file, refined_families)
    biome_results  = PARSE_BIOMES(query_results)
    domain_results = PARSE_DOMAINS(query_results, refined_families)

    db = INITIATE_SQLITE(params.db_schema_file, tables)
    APPEND_BLOBS_PARALLEL(db, cif_ch, seed_msa_ch, msa_ch, hmm_ch, rf_ch, \
        biome_results, domain_results)
}
