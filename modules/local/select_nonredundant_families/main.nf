process SELECT_NONREDUNDANT_FAMILIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(edges)
    tuple val(meta2), path(metadata, stageAs: "metadata/*")
    val(redundant_score_threshold)
    val(redundant_residue_threshold)

    output:
    tuple val(meta), path("${prefix}_redundant.txt"), emit: txt
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_nonredundant_families.py \\
        --edges ${edges} \\
        --metadata metadata/* \\
        --score_threshold ${redundant_score_threshold} \\
        --residue_threshold ${redundant_residue_threshold} \\
        --out ${prefix}_redundant.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_redundant.txt
    """
}
