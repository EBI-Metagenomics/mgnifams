process QUERY_MGNPROTEIN_DB {
    publishDir "${params.outDir}/post-processing", mode: "copy"
    label "venv"

    input:
    path config_file
    path family_proteins_file

    output:
    path "query_results"

    """
    python3 ${params.scriptDir}/post-processing/query_mgnprotein_db.py ${config_file} ${family_proteins_file}
    """
}

process PARSE_BIOMES {
    publishDir "${params.outDir}/post-processing", mode: "copy"
    label "venv"

    input:
    path query_results

    output:
    path "biome_results"

    """
    python3 ${params.scriptDir}/post-processing/parse_biomes.py ${query_results}
    """
}

process PARSE_DOMAINS {
    publishDir "${params.outDir}/post-processing", mode: "copy"
    label "venv"

    input:
    path query_results

    output:
    path "domain_results"

    """
    python3 ${params.scriptDir}/post-processing/parse_domains.py ${query_results}
    """
}
