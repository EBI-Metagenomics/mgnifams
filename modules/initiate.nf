process FILTER_UNANNOTATED_SLICES {
    label "venv"

    input:
    path sequence_chunk
    val min_sequence_length

    output:
    path "${sequence_chunk.baseName}.fa"

    script:
    """
    python3 ${params.scriptDir}/filter_unannotated_slices_fasta.py ${sequence_chunk} ${sequence_chunk.baseName}.fa ${min_sequence_length}
    """
}

process PUBLISH_INPUT_FASTA {
    input:
    path fasta

    output:
    path "${fasta}"

    script:
    """
    echo 'MGnifams input fasta published.'
    """
}