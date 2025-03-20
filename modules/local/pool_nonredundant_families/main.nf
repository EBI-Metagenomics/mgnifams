process POOL_NONREDUNDANT_FAMILIES {
    tag "$meta12.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

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
    tuple val(meta11), path(family_reps_ch , stageAs: "families_prepooled/family_reps/*")
    tuple val(meta12), path(logs_ch        , stageAs: "families_prepooled/logs/*")
    tuple val(meta13), path(redundant_family_ids)
    val(starting_id)

    output:
    tuple val(meta)  , path("families/seed_msa_sto/*")         , emit: seed_msa_sto
    tuple val(meta2) , path("families/msa_sto/*")              , emit: msa_sto
    tuple val(meta3) , path("families/hmm/*")                  , emit: hmm
    tuple val(meta4) , path("families/rf/*")                   , emit: rf
    tuple val(meta5) , path("families/domtblout/*")            , emit: domtblout
    tuple val(meta6) , path("families/refined_families.tsv")   , emit: tsv
    tuple val(meta7) , path("families/discarded_clusters.txt") , emit: discarded
    tuple val(meta8) , path("families/successful_clusters.txt"), emit: successful
    tuple val(meta9) , path("families/converged_families.txt") , emit: converged
    tuple val(meta10), path("families/family_metadata.csv")    , emit: metadata
    tuple val(meta11), path("families/family_reps.fasta")      , emit: family_reps
    tuple val(meta13), path("families/family_to_id.json")      , emit: id_mapping
    path "versions.yml"                                        , emit: versions

    script:
    """
    pool_nonredundant_families.py \\
        --input_dir families_prepooled \\
        --output_dir families \\
        --redundant ${redundant_family_ids} \\
        --iter ${starting_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
