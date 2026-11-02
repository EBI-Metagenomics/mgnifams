process FINALIZE_SQLITE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d41320ac1ca5f0a626982296d23ded50376a966b9e8240aa50dba2014a805bf5/data':
        'community.wave.seqera.io/library/sqlite:3.48.0--48957425ca78aa09' }"

    input:
    tuple val(meta), path(db, stageAs: 'input/db.sqlite3')
    path sql // e.g. assets/finalize_db.sql

    output:
    tuple val(meta), path("${prefix}.sqlite3"), emit: db
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Work on a copy (-resume safe); it is a throwaway build file, so no journal or fsync is needed
    cp -L "${db}" "${prefix}.sqlite3"
    {
        echo "PRAGMA journal_mode = OFF;"
        echo "PRAGMA synchronous = OFF;"
        echo "PRAGMA temp_store = MEMORY;"
        echo "PRAGMA cache_size = -2000000;"
        cat "${sql}"
    } | sqlite3 -bail "${prefix}.sqlite3" > /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -L "${db}" "${prefix}.sqlite3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}
