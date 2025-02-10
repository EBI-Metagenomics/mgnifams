#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV   } from "../modules/local/export_mgnifams_csv.nf"
include { QUERY_MGNPROTEIN_DB   } from "../modules/local/query_mgnprotein_db.nf"
include { PARSE_BIOMES          } from "../modules/local/parse_biomes.nf"
include { PARSE_DOMAINS         } from "../modules/local/parse_domains.nf"

workflow EXPORT_DATA {
    take:
    fam_metadata
    fam_converged
    refined_families
    pfam_hits
    foldseek_hits
    predict_scores
    
    main:
    tables = EXPORT_MGNIFAMS_CSV( fam_metadata, fam_converged, \
        refined_families, pfam_hits, foldseek_hits, predict_scores )
    query_results  = QUERY_MGNPROTEIN_DB( Channel.of( [ [id:"config"], params.db_config_file ] ), refined_families )
    
    // Chunking query results to run in parallel
    def counter = 0
    ch_query_results_batch = query_results.res
        .map{ meta, path ->
            path
        }
        .flatten()
        .collate( params.query_result_chunks )
        .map { batch -> 
            counter += 1
            [ [id:"batch_" + counter], batch ]
        }

    PARSE_BIOMES( ch_query_results_batch, query_results.biome_mapping.first() )
    PARSE_DOMAINS( ch_query_results_batch, query_results.pfam_mapping.first(), refined_families.first() )
}
