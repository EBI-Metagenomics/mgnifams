include { INSERT_SQLITE       } from '../../../modules/local/insert_sqlite/main.nf'
// include { QUERY_MGNPROTEIN_DB } from '../../../modules/local/query_mgnprotein_db/main'
// include { PARSE_BIOMES        } from '../../../modules/local/parse_biomes/main'
// include { PARSE_DOMAINS       } from '../../../modules/local/parse_domains/main'
// include { APPEND_SQLITE_BLOBS } from '../../../modules/local/append_sqlite_blobs/main.nf'

workflow UPDATE_DB {
    take:
    samplesheet
    
    main:
    ch_versions = Channel.empty()

    ch_queries = samplesheet
        .multiMap { meta, pipeline_results, db ->
            insert: [ meta, file("${pipeline_results.toUriString()}/table_data/*", checkIfExists: true), db ]
            update: [ meta, file("${pipeline_results.toUriString()}", checkIfExists: true) ]
        }
    //TODO more multiMaps for query/parse modules

    INSERT_SQLITE( ch_queries.insert )
    ch_versions = ch_versions.mix( INSERT_SQLITE.out.versions )

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

    // TODO append PARSE_BIOMES and PARSE_DOMAINS from above
    // APPEND_SQLITE_BLOBS( INSERT_SQLITE.out.db, ch_queries.update )
    // ch_versions = ch_versions.mix( APPEND_SQLITE_BLOBS.out.versions )

    emit:
    versions = ch_versions
}
