process POOL_PREREDUNDANT_FAMILIES_TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(tsv_ch, stageAs: "refined_families/*")

    output:
    tuple val(meta), path("${prefix}_preredundant_families.tsv"), emit: tsv
    path "versions.yml"                                         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pool_preredundant_families_tsv.py \\
        --input_dir refined_families \\
        --output_file ${prefix}_preredundant_families.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preredundant_families.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
