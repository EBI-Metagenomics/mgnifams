process HHSUITE_REFORMAT {
    tag "$meta.id"
    label 'process_single'
    // publishDir "${params.outdir}/a3m", mode: "copy"
    
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fa)
    val(informat)
    val(outformat)

    output:
    tuple val(meta), path("${meta.id}"), emit: fa
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    for file in ${fa}/*; do
        name=\$(basename \$file .${informat})

        ${moduleDir}/templates/reformat.pl \\
            $args \\
            ${informat} \\
            ${outformat} \\
            \$file \\
            ${prefix}/\$name.${outformat}
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}.${outformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
