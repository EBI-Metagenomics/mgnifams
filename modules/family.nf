process REFINE_FAMILIES {
    publishDir "${params.outdir}/families/", mode: "copy"
    label "venv"

    input:
    tuple val(meta), path(families_tsv)
    path(fasta)
    val(minimum_members)

    output:
    tuple val(meta), path("refined_families.tsv"), emit: tsv
    path("log.txt")                              , emit: log

    script:
    """
    python3 ${params.scriptDir}/refine_families.py ${families_tsv} ${fasta} ${minimum_members} refined_families.tsv > log.txt
    """
}

process FILTER_FAMILIES {
    input:
    path clust_tsv
    val threshold

    output:
    path "filtered_clusters.tsv", emit: filtered_clusters

    script:
    """
    ${params.scriptDir}/filter_families.sh ${clust_tsv} ${threshold} filtered_clusters.tsv
    """
}

process PARSE_FAMILIES {
    publishDir "${params.outdir}/families/", mode: "copy"
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
    python3 ${params.scriptDir}/parse_families.py ${fasta} ${clust_tsv}
    """
}
