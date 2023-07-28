process EXPORT_PROTEINS_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path fasta

    output:
    path "mgnifams_proteins.csv"

    """
    python3 ${baseDir}/bin/export_proteins_csv.py ${fasta} mgnifams_proteins.csv
    """
}

process EXPORT_CLUSTERING_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path cluster_file

    output:
    path "mgnifams_clusters.csv"

    """
    echo 'ProteinID,FamilyID' > mgnifams_clusters.csv && awk -v OFS="," '{print \$2,\$1}' ${cluster_file} >> mgnifams_clusters.csv
    """
}

process EXPORT_FAMILIES_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path rep_names
    val mode

    output:
    path "mgnifams_${mode}_families.csv"

    """
    python3 ${baseDir}/bin/export_families_csv.py ${rep_names} mgnifams_${mode}_families.csv ${mode}
    """
}