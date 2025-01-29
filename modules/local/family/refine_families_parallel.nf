process REFINE_FAMILIES_PARALLEL {
    maxForks 30
    conda "${moduleDir}/environment.yml"
    
    input:
    tuple val(meta), path(clusters_chunk)
    tuple val(meta2), path(mgnifams_fasta)

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
    path("versions.yml")         , topic: 'versions'

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
