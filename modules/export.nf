process EXPORT_PROTEINS_CSV { // TODO maybe also execute this afterwards, along with EXPORT_MGNIFAMS_CSV
    publishDir "${params.outdir}/tables", mode: "copy"
    label "venv"

    input:
    path fasta

    output:
    path "proteins.csv"

    """
    python3 ${params.scriptDir}/export_proteins_csv.py ${fasta} proteins.csv
    """
}

process EXPORT_MGNIFAMS_CSV {
    publishDir "${params.outdir}/tables", mode: "copy"
    label "venv"

    input:
    path out_dir

    output:
    path "mgnifams.csv"
    path "mgnifams_pfams.csv"
    path "mgnifams_folds.csv"

    """
    python3 ${params.scriptDir}/export_mgnifams_csvs.py ${out_dir}
    """
}