process EXPORT_MGNIFAMS_CSV {
    publishDir "${params.outDir}", mode: "copy"
    label "venv"

    input:
    path out_dir

    output:
    path "tables/mgnifam.csv"
    path "tables/mgnifam_proteins.csv"
    path "tables/mgnifam_pfams.csv"
    path "tables/mgnifam_folds.csv"

    """
    python3 ${params.scriptDir}/export_mgnifams_csvs.py ${out_dir}
    """
}
