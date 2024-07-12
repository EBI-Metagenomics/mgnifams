process POOL_FAMILY_RESULTS {
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
