process QUERY_MGNPROTEIN_DB {
    publishDir "${params.outDir}/post-processing", mode: "copy"
    label "venv"

    input:
    path config_file
    path family_proteins_file

    output:
    path "query_results"

    """
    python3 ${params.scriptDir}/query_mgnprotein_db.py ${config_file} ${family_proteins_file}
    """
}
