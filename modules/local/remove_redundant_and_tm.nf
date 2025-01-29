process REMOVE_REDUNDANT_AND_TM {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(hhblits_hits)
    tuple val(meta2), path(fam_rep_mapping)
    tuple val(meta3), path(tm_ids)
    tuple val(meta4), path(prob_ids)
    tuple val(meta5), path(refined_fam_proteins)
    tuple val(meta6), path(rep_fa)
    val redundant_threshold
    val similarity_threshold

    output:
    tuple val(meta), path("non_redundant_fam_ids.txt"), emit: non_redundant_fam_ids
    tuple val(meta), path("redundant_fam_ids.txt")    , emit: redundant_fam_ids
    tuple val(meta), path("similarity_edgelist.csv")  , emit: similarity_edgelist
    tuple val(meta), path("log.txt")                  , emit: log

    script:
    """
    remove_redundant_and_tm.py \\
        ${hhblits_hits} ${fam_rep_mapping} \\
        ${tm_ids} ${prob_ids} ${refined_fam_proteins} ${rep_fa} \\
        non_redundant_fam_ids.txt redundant_fam_ids.txt similarity_edgelist.csv log.txt \\
        ${redundant_threshold} ${similarity_threshold}
    """
}
