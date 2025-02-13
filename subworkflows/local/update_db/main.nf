include { INSERT_SQLITE       } from "../../../modules/local/insert_sqlite/main.nf"
include { APPEND_SQLITE_BLOBS } from "../../../modules/local/append_sqlite_blobs/main.nf"

workflow UPDATE_DB {
    take:
    samplesheet
    
    main:
    ch_versions = Channel.empty()

    ch_queries = samplesheet
        .multiMap { meta, pipeline_results, db ->
            insert: [ meta, file("${pipeline_results.toUriString()}/post-processing/tables/*", checkIfExists: true), db ]
            update: [ meta, file("${pipeline_results.toUriString()}", checkIfExists: true) ]
        }
    
    INSERT_SQLITE( ch_queries.insert )
    ch_versions = ch_versions.mix( INSERT_SQLITE.out.versions )
    APPEND_SQLITE_BLOBS( INSERT_SQLITE.out.db, ch_queries.update )
    ch_versions = ch_versions.mix( APPEND_SQLITE_BLOBS.out.versions )

    emit:
    versions = ch_versions
}
