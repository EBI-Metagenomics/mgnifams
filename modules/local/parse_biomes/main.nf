process PARSE_BIOMES {
    tag "$meta.id"
    
    input:
    tuple val(meta) , path(query_results, stageAs: "query_results/*")
    tuple val(meta2), path(biome_mapping)

    output:
    tuple val(meta), path("biome_results/*")

    script:
    """
    parse_biomes.py \\
        query_results \\
        ${biome_mapping}
    """
}
