include { BUILD_PARQUET_BIOME_QUERIES } from '../../../modules/local/build_parquet_biome_queries/main'
include { PARSE_BIOMES                } from '../../../modules/local/parse_biomes/main'
include { batchQueryResults           } from '../utils_nfcore_mgnifams_pipeline'

//
// Biome distributions of family members from the MGnify clusters parquet (cluster_rep_biomes): each member's own
// biomes, so members must be cluster representatives. Used by update_mgnifams; reusable by a biome-only rebuild.
//
workflow BIOMES_FROM_PARQUET {
    take:
    ch_family_chunks    // channel: [ meta, [ family<TAB>member TSVs ] ], whole families per chunk
    ch_clusters         // channel: clusters parquet (value)
    query_result_chunks // integer: families per PARSE_BIOMES task

    main:
    ch_versions = channel.empty()

    BUILD_PARQUET_BIOME_QUERIES( ch_family_chunks.combine( ch_clusters ) )
    ch_versions = ch_versions.mix( BUILD_PARQUET_BIOME_QUERIES.out.versions )

    ch_biome_batches = BUILD_PARQUET_BIOME_QUERIES.out.res
        .flatMap { meta, files -> batchQueryResults( meta, files, query_result_chunks ) }
    PARSE_BIOMES( ch_biome_batches, [ [id: 'lineages'], [] ] )
    ch_versions = ch_versions.mix( PARSE_BIOMES.out.versions )

    emit:
    biomes   = PARSE_BIOMES.out.res // channel: [ meta, biome_results/<family>.csv ]
    versions = ch_versions
}
