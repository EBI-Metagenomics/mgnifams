process EXPORT_MGNIFAMS_CSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta ), path(metadata)
    tuple val(meta2), path(scores)
    // tuple val(meta4), path(pfam_hits   , stageAs: "results/hh/*")
    // tuple val(meta5), path(foldseek_hit, stageAs: "results/structures/foldseek/*")
    // tuple val(meta6), path(scores      , stageAs: "results/structures/*")

    output:
    tuple val(meta), path("mgnifam.csv"), emit: csv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export_mgnifams_csvs.py \\
        --metadata ${metadata} \\
        --structure_scores ${scores} \\
        --mgnifam_out mgnifam.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """

    stub:
    """
    touch mgnifam.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """
}
