process APPEND_BLOBS_PARALLEL {
    tag "$meta.id"

    input:
    tuple val(meta) , path(db)
    tuple val(meta2), path(cif_ch     , stageAs: "output/structures/cif/*")
    tuple val(meta3), path(seed_msa_ch, stageAs: "output/families/*")
    tuple val(meta4), path(msa_ch     , stageAs: "output/families/*")
    tuple val(meta5), path(hmm_ch     , stageAs: "output/families/hmm/*")
    tuple val(meta6), path(rf_ch      , stageAs: "output/families/rf/*")
    tuple val(meta7), path(biome_ch   , stageAs: "output/post-processing/*")
    tuple val(meta8), path(domain_ch  , stageAs: "output/post-processing/*")
    
    output:
    tuple val(meta), path("${db}")

    script:
    """
    append_blobs_sqlite_parallel.py \\
        ${db} \\
        output \\
        families \\
        ${task.cpus}
    """
}
