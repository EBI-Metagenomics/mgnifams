process KEEP_UNANNOTATED {
    label "venv"

    input:
    path mgnify90_csv

    output:
    path "unannotated.fa"

    script:
    """
    python3 ${params.scriptDir}keep_unannotated_fasta.py ${mgnify90_csv} unannotated.fa
    """
}