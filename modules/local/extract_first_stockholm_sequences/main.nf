process EXTRACT_FIRST_STOCKHOLM_SEQUENCES {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_sto, stageAs: "msa_sto/*")

    output:
    tuple val(meta), path("family_reps.fasta")  , emit: fa
    tuple val(meta), path("problematic_ids.txt"), emit: prob_ids

    script:
    """
    extract_first_stockholm_sequences.py \\
        msa_sto \\
        family_reps.fasta \\
        problematic_ids.txt
    """
}
