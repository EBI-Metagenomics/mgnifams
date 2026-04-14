include { INIT_SQLITE    } from '../modules/local/init_sqlite/main'
include { IMPORT_QUERIES } from '../modules/local/import_queries/main.nf'

workflow INIT_DB {
    take:
    samplesheet
    
    main:
    ch_versions = Channel.empty()

    ch_queries = samplesheet
        .multiMap { meta, schema, pipeline_results ->
            schema: [ meta, file("${schema}", checkIfExists: true) ]
            pipeline_results: [ meta, file("${pipeline_results.toUriString()}/table_data/*", checkIfExists: true) ]
        }

    INIT_SQLITE( ch_queries.schema )
    ch_versions = ch_versions.mix( INIT_SQLITE.out.versions )

    IMPORT_QUERIES( ch_queries.pipeline_results, INIT_SQLITE.out.db )
    ch_versions = ch_versions.mix( IMPORT_QUERIES.out.versions )

    emit:
    versions = ch_versions
}
