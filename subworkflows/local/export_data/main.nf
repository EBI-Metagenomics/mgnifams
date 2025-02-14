#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV   } from "../../../modules/local/export_mgnifams_csv/main"
include { QUERY_MGNPROTEIN_DB   } from "../../../modules/local/query_mgnprotein_db/main"
include { PARSE_BIOMES          } from "../../../modules/local/parse_biomes/main"
include { PARSE_DOMAINS         } from "../../../modules/local/parse_domains/main"

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
    ch_query_results_batch = query_results.res
        .flatMap { _meta, files -> 
            files.collate( params.query_result_chunks )
                .withIndex()
                .collect{ flist, index -> tuple( [id: "batch_${index}" ], flist ) }
        }

    PARSE_BIOMES( ch_query_results_batch, query_results.biome_mapping.first() )
    PARSE_DOMAINS( ch_query_results_batch, query_results.pfam_mapping.first(), refined_families.first() )
}
