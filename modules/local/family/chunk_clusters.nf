process CHUNK_CLUSTERS {
    conda "${moduleDir}/environment.yml"

    input:
    path(clusters)
    path(checked_clusters)
    val(minimum_members)
    val(num_cluster_chunks)

    output:
    path("cluster_chunks/*")

    script:
    def checked_clusters = "${checked_clusters}" ?: 0
    """
    chunk_clusters.py ${clusters} ${checked_clusters} ${minimum_members} ${num_cluster_chunks}
    """
}
