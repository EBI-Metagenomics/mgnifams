process FOLDSEEK_EASYSEARCH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outDir}/foldseek", mode: "copy"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek:8.ef4e960--pl5321hb365157_0':
        'biocontainers/foldseek:8.ef4e960--pl5321hb365157_0' }"

    input:
    tuple val(meta)   , path(pdb)
    tuple val(meta_db), path(db)

    output:
    tuple val(meta), path("*.m8"), emit: aln
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    foldseek \\
        easy-search \\
        ${pdb} \\
        ${db}/${meta_db.id} \\
        ${meta_db.id}_${prefix}.m8 \\
        tmpFolder \\
        -e 0.001 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${meta_db.id}_${prefix}.m8

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}
