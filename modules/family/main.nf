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

process REFINE_FAMILIES_PARALLEL {
    conda "${moduleDir}/environment.yml"
    errorStrategy 'finish'
    
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
    path("logs/*")               , emit: logs

    script:
    """
    refine_families_parallel.py ${clusters_chunk} ${mgnifams_fasta} ${task.cpus}
    """
}

process MOVE_TO_DIR {
    label "venv"

    input:
    path(files)
    val(dir_name)

    output:
    path("${dir_name}")

    script:
    """
    mkdir -p ${dir_name}
    cp -r ${files} ${dir_name}
    """
}

process EXTRACT_FIRST_STOCKHOLM_SEQUENCES {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_sto, stageAs: "msa_sto/*")

    output:
    tuple val(meta), path("family_reps.fasta")

    script:
    """
    extract_first_stockholm_sequences.py msa_sto family_reps.fasta
    """
}

process EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_sto)

    output:
    path "family_reps.fasta"

    script:
    """
    extract_first_stockholm_sequences.py ${msa_sto} family_reps.fasta
    """
}

process MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(msa_a3m)

    output:
    path "fam_rep_mapping.csv"

    script:
    """
    map_first_a3m_sequences_to_family_id.py ${msa_a3m} fam_rep_mapping.csv
    """
}

process POOL_FAM_PROTEINS {
    label "venv"

    input:
    path tsv, stageAs: "refined_families/*"

    output:
    path "fam_proteins.tsv"

    script:
    """
    pool_fam_proteins.py refined_families fam_proteins.tsv
    """
}

process REMOVE_REDUNDANT_AND_TM {
    publishDir "${params.outdir}/redundancy/", mode: "copy"
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(hhblits_hits)
    path(fam_rep_mapping)
    tuple val(meta2), path(tm_ids)
    path(refined_fam_proteins)
    tuple val(meta3), path(rep_fa)

    output:
    path("non_redundant_fam_ids.txt"), emit: non_redundant_fam_ids
    path("redundant_fam_ids.txt")    , emit: redundant_fam_ids
    path("similarity_edgelist.csv")  , emit: similarity_edgelist
    path("log.txt")                  , emit: log

    script:
    """
    remove_redundant_and_tm.py \
        ${hhblits_hits} ${fam_rep_mapping} \
        ${tm_ids} ${refined_fam_proteins} ${rep_fa} \
        non_redundant_fam_ids.txt redundant_fam_ids.txt similarity_edgelist.csv log.txt
    """
}

process POOL_FAMILY_RESULTS {
    publishDir "${params.outdir}/", mode: "copy"
    conda "${moduleDir}/environment.yml"

    input:
    path seed_msa_sto_ch, stageAs: "families_prepooled/seed_msa_sto/*"
    path msa_sto_ch     , stageAs: "families_prepooled/msa_sto/*"
    path hmm_ch         , stageAs: "families_prepooled/hmm/*"
    path rf_ch          , stageAs: "families_prepooled/rf/*"
    path domtblout_ch   , stageAs: "families_prepooled/domtblout/*"
    path tsv_ch         , stageAs: "families_prepooled/refined_families/*"
    path discarded_ch   , stageAs: "families_prepooled/discarded_clusters/*"
    path successful_ch  , stageAs: "families_prepooled/successful_clusters/*"
    path converged_ch   , stageAs: "families_prepooled/converged_families/*"
    path metadata_ch    , stageAs: "families_prepooled/family_metadata/*"
    path logs_ch        , stageAs: "families_prepooled/logs/*"
    path non_redundant_family_ids
    path similarity_edgelist

    output:
    path("families/seed_msa_sto/*")         , emit: seed_msa_sto
    path("families/msa_sto/*")              , emit: msa_sto
    path("families/hmm/*")                  , emit: hmm
    path("families/rf/*")                   , emit: rf
    path("families/domtblout/*")            , emit: domtblout
    path("families/refined_families.tsv")   , emit: tsv
    path("families/discarded_clusters.txt") , emit: discarded
    path("families/successful_clusters.txt"), emit: successful
    path("families/converged_families.txt") , emit: converged
    path("families/family_metadata.csv")    , emit: metadata
    path("families/family_to_id.json")      , emit: id_mapping

    script:
    """
    pool_results.py families_prepooled families ${non_redundant_family_ids} ${similarity_edgelist}
    """
}
