process BUILD_PYFASTX_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyfastx:2.2.0--py39h0699b22_0':
        'biocontainers/pyfastx:2.2.0--py39h0699b22_0' }"
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}.fxi"), emit: index
    path("versions.yml")                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_pyfastx_index.py \\
        --fasta_file ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyfastx: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfastx'))")
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fxi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyfastx: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfastx'))")
    END_VERSIONS
    """
}
