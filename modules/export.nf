process EXPORT_PROTEINS_CSV {
    publishDir "${params.outdir}/tables", mode: "copy"
    label "venv"

    input:
    path fasta

    output:
    path "proteins.csv"

    """
    python3 ${params.scriptDir}/export/export_proteins_csv.py ${fasta} proteins.csv
    """
}

process EXPORT_MGNIFAMS_CSV {
    publishDir "${params.outdir}", mode: "copy"
    label "venv"

    input:
    path out_dir

    output:
    path "tables/mgnifam.csv"
    path "tables/mgnifam_proteins.csv"
    path "tables/mgnifam_pfams.csv"
    path "tables/mgnifam_folds.csv"

    """
    python3 ${params.scriptDir}/export/export_mgnifams_csvs.py ${out_dir}
    """
}