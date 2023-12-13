process HHSUITE_HHBLITS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hhsuite:3.3.0--py310pl5321h2a84d7f_5':
        'biocontainers/hhsuite:3.3.0--py310pl5321h2a84d7f_5' }"

    input:
    tuple val(meta), path(a3m_folder)
    path(hmm_db)
    val(db_name)

    output:
    tuple val(meta), path("${meta.id}"), emit: hhr
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    for a3m_file in ${a3m_folder}/*; do
        name=\$(basename \$a3m_file .a3m)
        hhblits \\
            $args \\
            -cpu $task.cpus \\
            -i \$a3m_file \\
            -d ${hmm_db}/${db_name} \\
            -o ${prefix}/\$name.hhr
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

    touch ${prefix}/${prefix}.hhr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
