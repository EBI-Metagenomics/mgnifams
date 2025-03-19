process REMOVE_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(domtbl)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(metadata, stageAs: "metadata/*")
    val(length_threshold)

    output:
    tuple val(meta), path("non_redundant.txt"), emit: txt
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    remove_redundant_fams.py \\
        --domtbl ${domtbl} \\
        --fasta ${fasta} \\
        --metadata metadata \\
        --length_threshold ${length_threshold} \\
        --out_file non_redundant.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """
}
