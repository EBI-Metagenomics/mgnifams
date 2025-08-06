process EXPORT_MGNIFAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta ), path(metadata)
    tuple val(meta2), path(scores)
    tuple val(meta3), path(composition)
    tuple val(meta4), path(tm_composition)

    output:
    tuple val(meta), path("mgnifam.csv"), emit: csv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def tm = tm_composition ? "${tm_composition}" : '""'
    """
    export_mgnifams.py \\
        --metadata ${metadata} \\
        --structure_scores ${scores} \\
        --composition ${composition} \\
        --tm_composition ${tm} \\
        --outfile mgnifam.csv
    
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
