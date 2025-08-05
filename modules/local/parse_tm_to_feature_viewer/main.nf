process PARSE_TM_TO_FEATURE_VIEWER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(preds) 

    output:
    tuple val(meta), path("${prefix}")         , emit: features
    tuple val(meta), path("tm_composition.csv"), emit: composition
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_tm_to_feature_viewer.py \\
        --input_file ${preds} \\
        --output_dir ${prefix} \\
        --csv_out tm_composition.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.json
    touch tm_composition.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
