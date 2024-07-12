process POOL_FAM_PROTEINS {
    label "venv"

    input:
    path tsv, stageAs: "refined_families/*"

    output:
    path "fam_proteins.tsv"

    script:
    """
    pool_fam_proteins.py refined_families fam_proteins.tsv
    """
}
