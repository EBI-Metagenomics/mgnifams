#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV   } from "../modules/local/export_mgnifams_csv.nf"
include { QUERY_MGNPROTEIN_DB   } from "../modules/local/query_mgnprotein_db.nf"
include { PARSE_BIOMES          } from "../modules/local/parse_biomes.nf"
include { PARSE_DOMAINS         } from "../modules/local/parse_domains.nf"
include { INITIATE_SQLITE       } from "../modules/local/initiate_sqlite.nf"
include { APPEND_BLOBS_PARALLEL } from "../modules/local/append_blobs_parallel.nf"

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
    fasta_input_mode
    
    main:
    tables = EXPORT_MGNIFAMS_CSV(fam_metadata, fam_converged, \
        refined_families, pfam_hits, foldseek_hits, predict_scores)
    if (!fasta_input_mode) {
        query_results  = QUERY_MGNPROTEIN_DB(params.db_config_file, refined_families)
        biome_results  = PARSE_BIOMES(query_results)
        domain_results = PARSE_DOMAINS(query_results, refined_families)

        db = INITIATE_SQLITE(params.db_schema_file, tables)
        APPEND_BLOBS_PARALLEL(db, cif_ch, seed_msa_ch, msa_ch, hmm_ch, rf_ch, \
            biome_results, domain_results)
    }
}
