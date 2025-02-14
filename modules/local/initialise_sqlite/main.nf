process INITIALISE_SQLITE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d41320ac1ca5f0a626982296d23ded50376a966b9e8240aa50dba2014a805bf5/data':
        'community.wave.seqera.io/library/sqlite:3.48.0--48957425ca78aa09' }"

    input:
    tuple val(meta), path(schema_file)
    
    output:
    tuple val(meta), path("${prefix}.sqlite3"), emit: db
    path "versions.yml"                       , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sqlite3 ${prefix}.sqlite3 < ${schema_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}
