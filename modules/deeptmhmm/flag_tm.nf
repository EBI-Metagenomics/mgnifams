process FLAG_TM {
    publishDir "${params.outDir}/tm", mode: "copy"
    tag "$meta.id"
    label 'venv'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gff3)
    val(fraction)

    output:
    tuple val(meta), path("tm_ids.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${params.scriptDir}/flag_transmembrane.py ${fraction} ${gff3} tm_ids.txt
    """
}
