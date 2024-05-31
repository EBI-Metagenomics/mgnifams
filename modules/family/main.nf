process CREATE_CLUSTERS_PKL {
    publishDir "${params.outDir}/families/", mode: "copy"
    
    conda "${moduleDir}/environment.yml"
    
    input:
    path(clusters_tsv)

    output:
    path("clusters_bookkeeping_df.pkl"), emit: pkl
    path("pkl_log.txt")                , emit: pkl_log
    path("refined_families.tsv")       , emit: refined_families
    path("discarded_clusters.txt")     , emit: discarded_clusters
    path("converged_families.txt")     , emit: converged_families
    path("family_metadata.csv")        , emit: family_metadata

    script:
    """
    python3 ${params.scriptDir}/family/create_clusters_bookkeeping_df.py ${clusters_tsv}
    touch refined_families.tsv
    touch discarded_clusters.txt
    touch converged_families.txt
    touch family_metadata.csv
    """
}

process REFINE_FAMILIES {
    publishDir "${params.outDir}/families/", mode: "copy"
    
    conda "${moduleDir}/environment.yml"
    
    input:
    path(clusters_pkl)
    path(families_tsv)
    path(mgnifams_fasta)
    path(discarded_clusters)
    path(converged_families)
    path(family_sizes)
    val(starting_num_sequences)
    val(minimum_members)
    val(iteration)

    output:
    path("seed_msa_sto/*")                , emit: seed_msa_sto
    path("msa_sto/*")                     , emit: msa_sto
    path("hmm/*")                         , emit: hmm
    path("rf/*")                          , emit: rf
    path("domtblout/*")                   , emit: domtblout
    path("updated_refined_families.tsv")  , emit: tsv
    path("updated_mgnifams_input.fa")     , emit: fa
    path("updated_discarded_clusters.txt"), emit: discarded
    path("updated_converged_families.txt"), emit: converged
    path("updated_family_metadata.csv")   , emit: metadata
    path("log.txt")                       , emit: log

    script:
    """
    python3 ${params.scriptDir}/family/refine_families.py ${clusters_pkl} ${families_tsv} ${mgnifams_fasta} ${discarded_clusters} ${converged_families} ${family_sizes} ${starting_num_sequences} ${minimum_members} ${iteration}
    """
}

process REFINE_FAMILIES_PARALLEL {
    publishDir "${params.outDir}/families/", mode: "copy"
    conda "${moduleDir}/environment.yml"
    
    input:
    path(clusters_chunk)
    path(mgnifams_fasta)

    output:
    path("seed_msa_sto/*")       , emit: seed_msa_sto
    path("msa_sto/*")            , emit: msa_sto
    path("hmm/*")                , emit: hmm
    path("rf/*")                 , emit: rf
    path("domtblout/*")          , emit: domtblout
    path("refined_families/*")   , emit: tsv
    path("discarded_clusters/*") , emit: discarded
    path("successful_clusters/*"), emit: successful
    path("converged_families/*") , emit: converged
    path("family_metadata/*")    , emit: metadata
    path("logs/*")               , emit: log

    script:
    """
    python3 ${params.scriptDir}/family/refine_families_parallel.py ${clusters_chunk} ${mgnifams_fasta} ${task.cpus}
    """
}

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
    """
    python3 ${params.scriptDir}/family/chunk_clusters.py ${clusters} ${checked_clusters} ${minimum_members} ${num_cluster_chunks}
    """
}

process EXTRACT_FIRST_STOCKHOLM_SEQUENCES {
    label "venv"

    input:
    tuple val(meta), path(msa_sto)

    output:
    path "family_reps.fasta"

    script:
    """
    python3 ${params.scriptDir}/extract_first_stockholm_sequences.py ${msa_sto} family_reps.fasta
    """
}
