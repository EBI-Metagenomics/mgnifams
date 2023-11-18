process FILTER_UNANNOTATED {
    label "venv"

    input:
    path mgnify90_csv

    output:
    path "mgnifams_input.fa"

    script:
    """
    python3 ${params.scriptDir}/filter_unannotated_fasta.py ${mgnify90_csv} mgnifams_input.fa
    """
}

process FILTER_UNANNOTATED_SLICES {
    label "venv"

    input:
    path mgnify90_csv
    val min_slice_length

    output:
    path "mgnifams_input.fa"

    script:
    """
    python3 ${params.scriptDir}/filter_unannotated_slices_fasta.py ${mgnify90_csv} mgnifams_input.fa ${min_slice_length}
    """
}
