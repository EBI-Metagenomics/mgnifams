process EXTRACT_LONG_FA {
    label "venv"

    input:
    path fasta , stageAs: "fasta/*"
    path scores, stageAs: "scores/*"

    output:
    path "family_reps.fasta"

    script:
    """
    python3 ${params.scriptDir}/extract_long_fa.py fasta scores family_reps.fasta
    """
}
