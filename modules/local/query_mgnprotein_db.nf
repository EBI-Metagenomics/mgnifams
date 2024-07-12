process QUERY_MGNPROTEIN_DB {
    label "venv"

    input:
    path config_file
    path family_proteins_file

    output:
    path "post-processing"

    """
    query_mgnprotein_db.py ${config_file} ${family_proteins_file}
    """
}
