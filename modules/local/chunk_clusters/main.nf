process CHUNK_CLUSTERS { // TODO repurpose
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::polars=1.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/19b2e332c6af2e0b5a319689bfc2236715ace874a279a486fd69ff3b8c72e1c2/data' :
        'community.wave.seqera.io/library/polars:1.21.0--8a6b364d2d19a70b' }"

    input:
    tuple val(meta), path(clusters)
    path(checked_clusters)
    val(minimum_members)
    val(num_cluster_chunks)

    output:
    tuple val(meta), path("cluster_chunks/*"), emit: fasta_chunks
    path "versions.yml"                      , emit: versions

    script:
    def omit_clusters = "${checked_clusters}" ?: 0
    """
    chunk_clusters.py \\
        ${clusters} \\
        ${omit_clusters} \\
        ${minimum_members} \\
        ${num_cluster_chunks}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('polars').version)")
    END_VERSIONS
    """
}
