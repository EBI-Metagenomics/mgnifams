process EXPORT_PROTEINS_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path fasta

    output:
    path "mgnifams_proteins.csv"

    """
    python3 ${baseDir}/bin/fasta_to_csv.py ${fasta} mgnifams_proteins.csv
    """
}