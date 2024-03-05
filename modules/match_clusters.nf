process MATCH_CLUSTERS {
    label "venv"

    input:
    tuple val(meta) , path(blastp_results)
    path clusters

    output:
    path "matched_cluster_reps.csv", emit: csv

    script:
    """
    python3 ${params.scriptDir}/match_clusters.py ${blastp_results} ${clusters}
    """
}