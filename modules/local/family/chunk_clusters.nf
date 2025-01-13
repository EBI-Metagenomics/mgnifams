process CHUNK_CLUSTERS {
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(clusters)
    path(checked_clusters)
    val(minimum_members)
    val(num_cluster_chunks)

    output:
    path("cluster_chunks/*")

    script:
    def omit_clusters = "${checked_clusters}" ?: 0
    """
    chunk_clusters.py ${clusters} ${omit_clusters} ${minimum_members} ${num_cluster_chunks}
    """
}
