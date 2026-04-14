process PARSE_DOMAINS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"
    
    input:
    tuple val(meta) , path(query_results, stageAs: "query_results/*")
    tuple val(meta2), path(pfam_mapping)
    tuple val(meta3), path(refined_families)

    output:
    tuple val(meta), path("domain_results/*"), emit: res
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_domains.py \\
        --query_results query_results \\
        --pfam_mapping ${pfam_mapping} \\
        --refined_families ${refined_families} \\
        --output_dir domain_results \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir domain_results
    touch domain_results/test.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """
}
