process EXTRACT_UNANNOTATED_SLICES {
    tag "$meta.id"

    conda "conda-forge::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb3700531c7ec639f59f084ab64c05e881d654dcf829db163539f2f0b095e09d/data' :
        'community.wave.seqera.io/library/biopython:1.84--3318633dad0031e7' }"

    input:
    tuple val(meta), path(sequence_chunk)
    val min_sequence_length

    output:
    tuple val(meta), path("${sequence_chunk.baseName}.fa"), emit: fa
    path "versions.yml"                                   , topic: 'versions'

    script:
    """
    extract_unannotated_slices.py \\
        ${sequence_chunk} \\
        ${sequence_chunk.baseName}.fa \\
        ${min_sequence_length}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
