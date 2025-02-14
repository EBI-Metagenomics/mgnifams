process MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_a3m)

    output:
    tuple val(meta), path("fam_rep_mapping.csv")

    script:
    """
    map_first_a3m_sequences_to_family_id.py \\
        ${msa_a3m} \\
        fam_rep_mapping.csv
    """
}
