process POOL_FAMILY_RESULTS {
    tag "$meta12.id"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta)  , path(seed_msa_sto_ch, stageAs: "families_prepooled/seed_msa_sto/*")
    tuple val(meta2) , path(msa_sto_ch     , stageAs: "families_prepooled/msa_sto/*")
    tuple val(meta3) , path(hmm_ch         , stageAs: "families_prepooled/hmm/*")
    tuple val(meta4) , path(rf_ch          , stageAs: "families_prepooled/rf/*")
    tuple val(meta5) , path(domtblout_ch   , stageAs: "families_prepooled/domtblout/*")
    tuple val(meta6) , path(tsv_ch         , stageAs: "families_prepooled/refined_families/*")
    tuple val(meta7) , path(discarded_ch   , stageAs: "families_prepooled/discarded_clusters/*")
    tuple val(meta8) , path(successful_ch  , stageAs: "families_prepooled/successful_clusters/*")
    tuple val(meta9) , path(converged_ch   , stageAs: "families_prepooled/converged_families/*")
    tuple val(meta10), path(metadata_ch    , stageAs: "families_prepooled/family_metadata/*")
    tuple val(meta11), path(logs_ch        , stageAs: "families_prepooled/logs/*")
    tuple val(meta12), path(non_redundant_family_ids)
    tuple val(meta13), path(similarity_edgelist)
    val(starting_id)

    output:
    tuple val(meta), path("families/seed_msa_sto/*")         , emit: seed_msa_sto
    tuple val(meta2), path("families/msa_sto/*")              , emit: msa_sto
    tuple val(meta3), path("families/hmm/*")                  , emit: hmm
    tuple val(meta4), path("families/rf/*")                   , emit: rf
    tuple val(meta5), path("families/domtblout/*")            , emit: domtblout
    tuple val(meta6), path("families/refined_families.tsv")   , emit: tsv
    tuple val(meta7), path("families/discarded_clusters.txt") , emit: discarded
    tuple val(meta8), path("families/successful_clusters.txt"), emit: successful
    tuple val(meta9), path("families/converged_families.txt") , emit: converged
    tuple val(meta10), path("families/family_metadata.csv")    , emit: metadata
    tuple val(meta13), path("families/similarity_edgelist.csv"), emit: similarity_edgelist
    tuple val(meta12), path("families/family_to_id.json")      , emit: id_mapping

    script:
    """
    pool_results.py \\
        families_prepooled \\
        families \\
        ${non_redundant_family_ids} \\
        ${similarity_edgelist} \\
        ${starting_id}
    """
}
