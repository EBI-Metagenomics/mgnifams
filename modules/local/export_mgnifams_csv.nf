process EXPORT_MGNIFAMS_CSV {
    tag "$meta.id"

    input:
    tuple val(meta),  path(metadata    , stageAs: "results/families/*")
    tuple val(meta2), path(converged   , stageAs: "results/families/*")
    tuple val(meta3), path(tsv         , stageAs: "results/families/*")
    tuple val(meta4), path(pfam_hits   , stageAs: "results/hh/*")
    tuple val(meta5), path(foldseek_hit, stageAs: "results/structures/foldseek/*")
    tuple val(meta6), path(scores      , stageAs: "results/structures/*")

    output:
    tuple val(meta), path("tables")

    script:
    """
    export_mgnifams_csvs.py \\
        results \\
        families \\
        tables
    """
}
