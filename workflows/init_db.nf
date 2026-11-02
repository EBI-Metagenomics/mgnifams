include { INIT_SQLITE     } from '../modules/local/init_sqlite/main'
include { IMPORT_QUERIES  } from '../modules/local/import_queries/main.nf'
include { FINALIZE_SQLITE } from '../modules/local/finalize_sqlite/main'

workflow INIT_DB {
    take:
    samplesheet

    main:
    ch_versions = channel.empty()

    ch_queries = samplesheet
        .multiMap { meta, schema, pipeline_results ->
            schema: [ meta, file("${schema}", checkIfExists: true) ]
            pipeline_results: [ meta, file("${pipeline_results.toUriString()}/table_data/*", checkIfExists: true) ]
        }

    INIT_SQLITE( ch_queries.schema )
    ch_versions = ch_versions.mix( INIT_SQLITE.out.versions )

    IMPORT_QUERIES( ch_queries.pipeline_results, INIT_SQLITE.out.db )
    ch_versions = ch_versions.mix( IMPORT_QUERIES.out.versions )

    // has_* flags, indexes, ANALYZE
    FINALIZE_SQLITE( IMPORT_QUERIES.out.db, file("${projectDir}/assets/finalize_db.sql", checkIfExists: true) )
    ch_versions = ch_versions.mix( FINALIZE_SQLITE.out.versions )
}
