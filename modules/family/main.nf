process REFINE_FAMILIES {
    publishDir "${params.outdir}/families/", mode: "copy"
    
    conda "${moduleDir}/environment.yml"
    
    input:
    path(families_tsv)
    path(fasta)
    val(minimum_members)
    path(clusters_bookkeeping_df_pkl)
    // path(updated_mgnifams_input_pkl)

    output:
    path("refined_families.tsv"), emit: tsv
    path("seed_msa/*")          , emit: seed_msa
    path("msa/*")               , emit: msa
    path("hmm/*")               , emit: hmm
    path("domtblout/*")         , emit: domtblout
    path("log.txt")             , emit: log

    script:
    """
    python3 ${params.scriptDir}/refine_families.py ${families_tsv} ${fasta} ${minimum_members} refined_families.tsv ${clusters_bookkeeping_df_pkl}
    """
    // python3 ${params.scriptDir}/refine_families.py ${families_tsv} ${fasta} ${minimum_members} refined_families.tsv
    // python3 ${params.scriptDir}/refine_families.py ${families_tsv} ${fasta} ${minimum_members} refined_families.tsv ${clusters_bookkeeping_df_pkl} ${updated_mgnifams_input_pkl}
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
