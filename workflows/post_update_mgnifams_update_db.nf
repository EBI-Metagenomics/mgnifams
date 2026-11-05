/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    POST_UPDATE_MGNIFAMS_UPDATE_DB: delta sqlite DB from an update_mgnifams or annotate_not_updated outdir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Only the successful families (update_families/family_ids.fasta), with the tables and
    blobs the update run published, ready to be merged (assets/sql/merge_update_delta.sql).
    UPDATE_SQLITE_BLOBS_STAGED fails instead of publishing a partial DB.
----------------------------------------------------------------------------------------
*/

include { INIT_SQLITE                } from '../modules/local/init_sqlite/main'
include { IMPORT_QUERIES             } from '../modules/local/import_queries/main'
include { UPDATE_SQLITE_BLOBS_STAGED } from '../modules/local/update_sqlite_blobs_staged/main'

workflow POST_UPDATE_MGNIFAMS_UPDATE_DB {

    take:
    samplesheet // channel: [ meta, update_mgnifams outdir, biome_results folder or [] ]

    main:
    ch_versions = channel.empty()

    // An annotate_not_updated outdir (annotation_only,true) has no structures, domains or biomes: prod keeps them
    ch_input = samplesheet
        .map { meta, results, biomes -> [ meta, results.toUriString(), biomes ? biomes.toUriString() : '' ] }
        .map { meta, r, b -> [ meta, r, b, file("${r}/update_families/update_info.csv").exists() && file("${r}/update_families/update_info.csv").readLines().contains('annotation_only,true') ] }
        .map { meta, r, b, annotation_only ->
            if (b && annotation_only) {
                error("post_update_mgnifams_update_db: biome_results is for update deltas; an annotation-only delta keeps prod's biomes")
            }
            [ meta, r, b, annotation_only ]
        }
        .multiMap { meta, r, b, annotation_only ->
            schema : [ meta, file("${projectDir}/assets/sql/db_schema.sqlite", checkIfExists: true) ]
            tables : [ meta, files("${r}/table_data/*.csv", checkIfExists: true) ]
            cifs   : annotation_only ? [] : files("${r}/structures/esmfold/cif/*.cif", checkIfExists: true)
            s4pred : files("${r}/annotation/reps/s4pred/json/*.json", checkIfExists: true)
            tm     : files("${r}/annotation/reps/deeptmhmm/json/*.json")
            biomes : annotation_only ? [] : files("${b ?: "${r}/biome_results"}/*")
            domains: annotation_only ? [] : files("${r}/domain_results/*", checkIfExists: true)
            ids    : file("${r}/update_families/family_ids.fasta", checkIfExists: true)
            // Written by update_mgnifams; also tells an update outdir apart from a run_mgnifams_pipeline one
            info   : [ file("${r}/update_families/update_info.csv", checkIfExists: true), b as boolean ]
        }

    INIT_SQLITE( ch_input.schema )
    ch_versions = ch_versions.mix( INIT_SQLITE.out.versions )

    IMPORT_QUERIES( ch_input.tables, INIT_SQLITE.out.db )
    ch_versions = ch_versions.mix( IMPORT_QUERIES.out.versions )

    ch_ids = ch_input.ids
        .map { fa -> fa.readLines().findAll { line -> line.startsWith('>') }.collect { line -> line.substring(1) }.sort().join('\n') }
        .collectFile( name: 'successful_ids.txt', newLine: true )

    // tm_computed / biome_computed make those blobs required; the merge reads them back from update_info
    ch_info = ch_input.info
        .map { csv, biomes_override ->
            def info = csv.splitCsv( header: true ).collectEntries { row -> [ row.key, row.value ] }
            info + [
                tm_computed   : info.tm_computed == 'true',
                // Rebuilt biomes (biome_results column) are computed whatever the update run recorded
                biome_computed: info.biome_computed == 'true' || biomes_override,
                child_tables  : 'pfams,funfams,folds'
            ] + ( info.annotation_only == 'true' ? [ annotation_only: true ] : [:] )
        }

    UPDATE_SQLITE_BLOBS_STAGED(
        IMPORT_QUERIES.out.db,
        ch_input.cifs,
        ch_input.s4pred,
        ch_input.tm,
        ch_input.biomes,
        ch_input.domains,
        ch_ids,
        ch_info
    )
    ch_versions = ch_versions.mix( UPDATE_SQLITE_BLOBS_STAGED.out.versions )
}
