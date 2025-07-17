include { EXPORT_MGNIFAMS   } from '../../../modules/local/export_mgnifams/main'
// include { EXPORT_FUNFAMS    } from '../../../modules/local/export_funfams/main'
// include { EXPORT_PFAMS      } from '../../../modules/local/export_pfams/main'
// include { EXPORT_FOLDS      } from '../../../modules/local/export_folds/main'

// include { QUERY_MGNPROTEIN_DB } from '../../../modules/local/query_mgnprotein_db/main'
// include { PARSE_BIOMES        } from '../../../modules/local/parse_biomes/main'
// include { PARSE_DOMAINS       } from '../../../modules/local/parse_domains/main'

workflow EXPORT_DATA {
    take:
    family_metadata
    predicted_scores
    composition
    funfam_domains
    pfam_hits
    foldseek_hits

    main:
    ch_versions = Channel.empty()

    EXPORT_MGNIFAMS( family_metadata, predicted_scores, composition )
    ch_versions = ch_versions.mix( EXPORT_MGNIFAMS.out.versions )

    // EXPORT_FUNFAMS( funfam_domains )
    // ch_versions = ch_versions.mix( EXPORT_FUNFAMS.out.versions )

    // EXPORT_PFAMS( pfam_hits )
    // ch_versions = ch_versions.mix( EXPORT_PFAMS.out.versions )

    // EXPORT_FOLDS( foldseek_hits )
    // ch_versions = ch_versions.mix( EXPORT_FOLDS.out.versions )


    // TODO move to own workflow
    // query_results  = QUERY_MGNPROTEIN_DB( Channel.of( [ [id:"config"], params.db_config_file ] ), refined_families )
    // // TODO ch_versions = ch_versions.mix( QUERY_MGNPROTEIN_DB.out.versions )
    
    // // Chunking query results to run in parallel
    // ch_query_results_batch = query_results.res
    //     .flatMap { _meta, files -> 
    //         files.collate( params.query_result_chunks )
    //             .withIndex()
    //             .collect{ flist, index -> tuple( [id: "batch_${index}" ], flist ) }
    //     }

    // PARSE_BIOMES( ch_query_results_batch, query_results.biome_mapping.first() )
    // // TODO ch_versions = ch_versions.mix( PARSE_BIOMES.out.versions )

    // PARSE_DOMAINS( ch_query_results_batch, query_results.pfam_mapping.first(), refined_families.first() )
    // // TODO ch_versions = ch_versions.mix( PARSE_DOMAINS.out.versions )

    emit:
    versions = ch_versions
}
