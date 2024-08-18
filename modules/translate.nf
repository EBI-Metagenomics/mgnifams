process TRANSLATE_MSA_MGYPS {
    input:
    tuple val(meta), path(msa_fas)

    output:
    tuple val(meta), path("${meta.id}"), emit: fa

    script:
    """
    python3 ${params.scriptDir}/translate_msa_mgyps.py ${msa_fas} ${meta.id}
    """
}
