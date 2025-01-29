process EXTRACT_LONG_FA {
    label "venv"

    input:
    tuple val(meta) , path(fasta , stageAs: "fasta/*")
    tuple val(meta2), path(scores, stageAs: "scores/*")

    output:
    tuple val(meta), path("family_reps.fasta")

    script:
    """
    extract_long_fa.py \\
        fasta \\
        scores \\
        family_reps.fasta
    """
}
