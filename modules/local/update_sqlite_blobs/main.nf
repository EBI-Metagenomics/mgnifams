process UPDATE_SQLITE_BLOBS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta) , path(db), path(pipeline_results)
    tuple val(meta2), path(biome_results, stageAs: "biome_results/*")
    tuple val(meta3), path(domain_results, stageAs: "domain_results/*")

    output:
    tuple val(meta), path(db), emit: db
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    update_sqlite_blobs.py \\
        --db ${db} \\
        --pipeline_results ${pipeline_results} \\
        --biome_results biome_results \\
        --domain_results domain_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch test.sqlite3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
