process QUERY_MGNPROTEIN_DB {
    container "quay.io/microbiome-informatics/mgnifams:2.0.0"

    input:
    tuple val(meta) , path(config_file)
    tuple val(meta2), path(family_proteins_file)

    output:
    tuple val(meta), path("post-processing")

    script:
    """
    query_mgnprotein_db.py \\
        ${config_file} \\
        ${family_proteins_file}
    """
}
