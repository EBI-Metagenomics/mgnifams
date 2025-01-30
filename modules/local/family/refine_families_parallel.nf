process REFINE_FAMILIES_PARALLEL {
    maxForks 30
    conda "${moduleDir}/environment.yml"
    
    input:
    tuple val(meta), path(clusters_chunk)
    tuple val(meta2), path(mgnifams_fasta)

    output:
    tuple val(meta), path("seed_msa_sto/*")       , emit: seed_msa_sto
    tuple val(meta), path("msa_sto/*")            , emit: msa_sto
    tuple val(meta), path("hmm/*")                , emit: hmm
    tuple val(meta), path("rf/*")                 , emit: rf
    tuple val(meta), path("domtblout/*")          , emit: domtblout
    tuple val(meta), path("refined_families/*")   , emit: tsv
    tuple val(meta), path("discarded_clusters/*") , emit: discarded
    tuple val(meta), path("successful_clusters/*"), emit: successful
    tuple val(meta), path("converged_families/*") , emit: converged
    tuple val(meta), path("family_metadata/*")    , emit: metadata
    tuple val(meta), path("logs/*")               , emit: logs
    path("versions.yml")                          , topic: 'versions'

    script:
    """
    refine_families_parallel.py \\
        ${clusters_chunk} \\
        ${mgnifams_fasta} \\
        ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        hmmer: \$(echo \$(hmmbuild -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
        easel: \$(esl-reformat -h | grep -o '^# Easel [0-9.]*' | sed 's/^# Easel *//')
    END_VERSIONS
    """
}
