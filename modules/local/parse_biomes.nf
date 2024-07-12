process PARSE_BIOMES {
    label "venv"

    input:
    path query_results

    output:
    path "biome_results"

    """
    parse_biomes.py ${query_results}
    """
}
