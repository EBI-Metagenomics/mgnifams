process EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_sto)

    output:
    path "family_reps.fasta"

    script:
    """
    extract_first_stockholm_sequences.py ${msa_sto} family_reps.fasta
    """
}
