process EXTRACT_UNANNOTATED_SLICES {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(sequence_chunk)
    val min_sequence_length

    output:
    tuple val(meta), path("${prefix}.fa"), emit: fa
    path "versions.yml"                  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_unannotated_slices.py \\
        --input_file ${sequence_chunk} \\
        --output_file ${prefix}.fa \\
        --min_sequence_length ${min_sequence_length}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
