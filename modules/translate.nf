process TRANSLATE_MSA_MGYPS {
    label "venv"

    input:
    tuple val(meta), path(msa_fas)

    output:
    tuple val(meta), path("${msa_fas}"), emit: fa

    script:
    """
    python3 ${params.scriptDir}/translate_msa_mgyps.py ${msa_fas}
    """
}
