process TRUNCATE_A3M {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(a3m, stageAs: 'input/*')
    val max_seqs

    output:
    tuple val(meta), path("${meta.id}.a3m.gz"), emit: a3m
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Keeps the first max_seqs records: the full MSA's representative plus its top-ranked recruits
    """
    set -euo pipefail

    gzip -cdf "${a3m}" | awk -v max=${max_seqs} '/^>/ { n++ } n <= max' | gzip -n > "${meta.id}.a3m.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1 | cut -d, -f1)
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip -n > "${meta.id}.a3m.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1 | cut -d, -f1)
    END_VERSIONS
    """
}
