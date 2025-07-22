include { QUERY_MGNPROTEIN_DB } from '../../../modules/local/query_mgnprotein_db/main'
include { PARSE_BIOMES        } from '../../../modules/local/parse_biomes/main'
include { PARSE_DOMAINS       } from '../../../modules/local/parse_domains/main'
// include { APPEND_SQLITE_BLOBS } from '../../../modules/local/append_sqlite_blobs/main.nf'

workflow UPDATE_DB {
    take:
    samplesheet
    query_result_chunks
    
    main:
    ch_versions = Channel.empty()

    ch_queries = samplesheet
        .multiMap { meta, pipeline_results, db, secrets ->
            query_mgnprotein: [ meta, secrets, file("${pipeline_results.toUriString()}/generate_families/families/refined_families.tsv", checkIfExists: true) ]
            refined_families: [ meta, file("${pipeline_results.toUriString()}/generate_families/families/refined_families.tsv", checkIfExists: true) ]
            update: [ meta, db, file("${pipeline_results.toUriString()}", checkIfExists: true) ]
        }

    QUERY_MGNPROTEIN_DB( ch_queries.query_mgnprotein )
    ch_versions = ch_versions.mix( QUERY_MGNPROTEIN_DB.out.versions )
    
    // Chunking query results to run in parallel
    ch_query_results_batch = QUERY_MGNPROTEIN_DB.out.res
        .flatMap { _meta, files -> 
            files.collate( query_result_chunks )
                .withIndex()
                .collect{ flist, index -> tuple( [id: "batch_${index}" ], flist ) }
        }

    PARSE_BIOMES( ch_query_results_batch, QUERY_MGNPROTEIN_DB.out.biome_mapping.first() )
    ch_versions = ch_versions.mix( PARSE_BIOMES.out.versions )

    PARSE_DOMAINS( ch_query_results_batch, QUERY_MGNPROTEIN_DB.out.pfam_mapping.first(), ch_queries.refined_families.first() )
    ch_versions = ch_versions.mix( PARSE_DOMAINS.out.versions )

    // TODO append PARSE_BIOMES and PARSE_DOMAINS from above
    // APPEND_SQLITE_BLOBS( ch_queries.update, PARSE_BIOMES.out., PARSE_DOMAINS.out., )
    // ch_versions = ch_versions.mix( APPEND_SQLITE_BLOBS.out.versions )

    emit:
    versions = ch_versions
}
