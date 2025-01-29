process PARSE_BIOMES {
    input:
    tuple val(meta), path(query_results)

    output:
    tuple val(meta), path("biome_results")

    script:
    """
    parse_biomes.py ${query_results}
    """
}
