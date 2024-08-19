process EXPORT_MGNIFAMS_CSV {
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
    export_mgnifams_csvs.py results families tables
    """
}
