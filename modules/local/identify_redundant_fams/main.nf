process IDENTIFY_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(domtbl)
    tuple val(meta2), path(metadata, stageAs: "metadata/*")
    tuple val(meta3), path(edgelist)
    val(length_threshold)
    val(redundant_score_threshold)
    val(similarity_score_threshold)

    output:
    tuple val(meta), path("${prefix}_redundant.txt")   , emit: txt
    tuple val(meta), path("${prefix}_similarities.csv"), emit: csv
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    identify_redundant_fams.py \\
        --domtbl ${domtbl} \\
        --metadata metadata \\
        --edgelist ${edgelist} \\
        --length_threshold ${length_threshold} \\
        --redundant_score_threshold ${redundant_score_threshold} \\
        --similarity_score_threshold ${similarity_score_threshold} \\
        --out_file ${prefix}_redundant.txt \\
        --similarity_csv ${prefix}_similarities.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_redundant.txt
    touch ${prefix}_similarities.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """
}
