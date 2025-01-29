process PARSE_BIOMES {
    input:
    path query_results

    output:
    path "biome_results"

    script:
    """
    parse_biomes.py ${query_results}
    """
}
