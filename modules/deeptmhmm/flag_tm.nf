process FLAG_TM {
    publishDir "${params.outDir}/redundancy/tm", mode: "copy"
    tag "$meta.id"
    label 'venv'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybiolib:1.1.2025--pyhdfd78af_0':
        'biocontainers/pybiolib:1.1.2025--pyhdfd78af_0' }"

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
