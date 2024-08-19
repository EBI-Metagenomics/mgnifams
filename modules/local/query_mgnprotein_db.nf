process QUERY_MGNPROTEIN_DB {
    container "quay.io/microbiome-informatics/mgnifams:2.0.0"

    input:
    path config_file
    path family_proteins_file

    output:
    path "post-processing"

    """
    query_mgnprotein_db.py ${config_file} ${family_proteins_file}
    """
}
