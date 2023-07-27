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

process EXPORT_FAMILIES_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path rep_names

    output:
    path "mgnifams_families.csv"

    """
    python3 ${baseDir}/bin/export_families_csv.py ${rep_names} mgnifams_families.csv
    """
}