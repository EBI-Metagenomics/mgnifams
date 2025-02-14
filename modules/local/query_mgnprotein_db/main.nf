process QUERY_MGNPROTEIN_DB {
    tag "$meta.id"

    container "quay.io/microbiome-informatics/mgnifams:2.0.0"

    input:
    tuple val(meta) , path(config_file)
    tuple val(meta2), path(family_proteins_file)

    output:
    tuple val(meta2), path("query_results/*")  , emit: res
    tuple val(meta2), path("biome_mapping.tsv"), emit: biome_mapping
    tuple val(meta2), path("pfam_mapping.tsv") , emit: pfam_mapping

    script:
    """
    query_mgnprotein_db.py \\
        ${config_file} \\
        ${family_proteins_file}
    """
}
