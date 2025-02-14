process TRANSLATE_MSA_MGYPS {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_fas)

    output:
    tuple val(meta), path("${meta.id}"), emit: fa

    script:
    """
    translate_msa_mgyps.py \\
        ${msa_fas} \\
        ${meta.id}
    """
}
