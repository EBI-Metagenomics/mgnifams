process EXPORT_MGNIFAMS_CSV {
    publishDir "${params.outDir}", mode: "copy"
    label "venv"

    input:
    path metadata                  , stageAs: "results/families/*"
    path converged                 , stageAs: "results/families/*"
    path tsv                       , stageAs: "results/families/*"
    tuple val(meta), path(pfam_hits, stageAs: "results/hh/*")
    path foldseek_hits             , stageAs: "results/structures/foldseek/*"
    path scores                    , stageAs: "results/structures/*"

    output:
    path "tables"

    """
    python3 ${params.scriptDir}/export_mgnifams_csvs.py results families tables
    """
}

process EXPORT_MGNIFAMS_CSV_FROM_FOLDER {
    publishDir "${params.outDir}", mode: "copy"
    label "venv"

    input:
    path out_dir

    output:
    path "tables"

    """
    python3 ${params.scriptDir}/export_mgnifams_csvs.py ${out_dir} families tables
    """
}
