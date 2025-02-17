process BZ2_UNCOMPRESS {
    tag "$meta.id"
    
    conda "conda-forge::bzip2=1.0.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e9d99a4ff7ab85c6e41671ef54982250c8cb847da90ba7221076d4ea7e23d02/data' :
        'community.wave.seqera.io/library/bzip2:1.0.8--6aa6fbbde430998d' }"

    input:
    tuple val(meta), path(bz2_file)

    output:
    tuple val(meta), path("${bz2_file.baseName}"), emit: file
    path "versions.yml"                          , emit: versions

    script:
    """
    bzip2 -d < ${bz2_file} > ${bz2_file.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bzip2: \$(echo \$(bzip2 --help 2>&1 | head -n 1 | sed 's/.*Version \\([0-9\\.]*\\).*/\\1/') )
    END_VERSIONS
    """
}
