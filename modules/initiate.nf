process FILTER_UNANNOTATED_SLICES {
    label "venv"

    input:
    path mgnify90_csv
    val min_slice_length

    output:
    path "${mgnify90_csv.baseName}.fa"

    script:
    """
    python3 ${params.scriptDir}/filter_unannotated_slices_fasta.py ${mgnify90_csv} ${mgnify90_csv.baseName}.fa ${min_slice_length}
    """
}

process PUBLISH_INPUT_FASTA {
    publishDir "${params.outDir}/input", mode: "copy"

    input:
    path fasta

    output:
    path "${fasta}"

    script:
    """
    echo MGnifams input fasta published.
    """
}