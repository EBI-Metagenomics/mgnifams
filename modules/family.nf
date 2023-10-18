process KEEP_FAMILIES {
    input:
    path clust_tsv
    val threshold

    output:
    path "filtered_clusters.tsv", emit: filtered_clusters

    script:
    """
    ${params.scriptDir}filter_families.sh ${clust_tsv} ${threshold} filtered_clusters.tsv
    """
}

process PARSE_FAMILIES {
    publishDir "${params.outdir}families/", mode: "copy"
    label "venv"

    input:
    path fasta
    path clust_tsv

    output:
    path "rep_names.txt", emit: reps_ids
    path "reps.fa"      , emit: reps_fasta
    path "families/*"   , emit: families_folder

    script:
    """
    python3 ${params.scriptDir}parse_families.py ${fasta} ${clust_tsv}
    """
}
