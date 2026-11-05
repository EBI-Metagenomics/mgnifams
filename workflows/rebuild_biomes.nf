/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REBUILD_BIOMES: biome distributions of a finished update's families, from the clusters parquet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The same BIOMES_FROM_PARQUET as update_mgnifams, on the update's published member list
    (update_families/refined_families.tsv). Publishes <outdir>/biome_results/<id>.csv; the
    biome_results samplesheet column of post_update_mgnifams_update_db puts them in the delta.
----------------------------------------------------------------------------------------
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SPLIT_FAMILY_MEMBERS   } from '../modules/local/split_family_members/main'
include { BIOMES_FROM_PARQUET    } from '../subworkflows/local/biomes_from_parquet'

workflow REBUILD_BIOMES {

    take:
    ch_samplesheet      // channel: [ meta, refined_families.tsv, clusters parquet ]
    family_query_chunks // integer: about this many lookup tasks, whole families each
    query_result_chunks // integer: families per PARSE_BIOMES task
    outdir

    main:
    ch_versions = channel.empty()

    SPLIT_FAMILY_MEMBERS( ch_samplesheet.map { meta, families, _clusters -> [ meta, families ] }, family_query_chunks )
    ch_versions = ch_versions.mix( SPLIT_FAMILY_MEMBERS.out.versions )

    BIOMES_FROM_PARQUET(
        SPLIT_FAMILY_MEMBERS.out.chunks.flatMap { _meta, chunks -> [chunks].flatten().collect { chunk -> [ [id: chunk.baseName], [ chunk ] ] } },
        ch_samplesheet.map { _meta, _families, clusters -> clusters }.first(),
        query_result_chunks
    )
    ch_versions = ch_versions.mix( BIOMES_FROM_PARQUET.out.versions )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML( ch_versions )
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )
}
