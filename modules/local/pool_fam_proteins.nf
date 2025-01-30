process POOL_FAM_PROTEINS {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(tsv, stageAs: "refined_families/*")

    output:
    tuple val(meta), path("fam_proteins.tsv")

    script:
    """
    pool_fam_proteins.py \\
        refined_families \\
        fam_proteins.tsv
    """
}
