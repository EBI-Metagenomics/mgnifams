process APPEND_BLOBS_PARALLEL {
    label "venv"

    input:
    path db
    path cif_ch   , stageAs: "output/structures/cif/*"
    tuple val(meta2), path(seed_msa_ch, stageAs: "output/families/*")
    tuple val(meta3), path(msa_ch     , stageAs: "output/families/*")
    path hmm_ch   , stageAs: "output/families/hmm/*"
    path rf_ch    , stageAs: "output/families/rf/*"
    path biome_ch , stageAs: "output/post-processing/*"
    path domain_ch, stageAs: "output/post-processing/*"
    
    output:
    path "${db}"

    """
    append_blobs_sqlite_parallel.py ${db} output families ${task.cpus}
    """
}
