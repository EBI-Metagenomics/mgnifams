process EXPORT_MGNIFAMS_CSV {
    publishDir "${params.outDir}", mode: "copy"
    label "venv"

    input:
    path out_dir

    output:
    path "tables"

    """
    python3 ${params.scriptDir}/export_mgnifams_csvs.py ${out_dir} families
    """
}
