process PARSE_BIOMES {
    input:
    path query_results

    output:
    path "biome_results"

    """
    parse_biomes.py ${query_results}
    """
}
