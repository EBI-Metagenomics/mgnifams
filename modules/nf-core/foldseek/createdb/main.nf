process FOLDSEEK_CREATEDB {
    tag "$meta.id"
    // label 'process_single'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0':
    //     'biocontainers/foldseek:9.427df8a--pl5321hb365157_0' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${meta.id}"), emit: db
    path "versions.yml"                , topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    foldseek \\
        createdb \\
        ${pdb} \\
        ${prefix}/${prefix} \\
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
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}
    touch ${prefix}/${prefix}_ca
    touch ${prefix}/${prefix}_ca.dbtype
    touch ${prefix}/${prefix}_ca.index
    touch ${prefix}/${prefix}_h
    touch ${prefix}/${prefix}_h.dbtype
    touch ${prefix}/${prefix}_h.index
    touch ${prefix}/${prefix}_ss
    touch ${prefix}/${prefix}_ss.dbtype
    touch ${prefix}/${prefix}_ss.index
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index
    touch ${prefix}/${prefix}.lookup
    touch ${prefix}/${prefix}.source

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek --help | grep Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}
