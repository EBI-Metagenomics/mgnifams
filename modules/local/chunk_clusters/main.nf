process CHUNK_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta) , path(reps)
    tuple val(meta2), path(clusters)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    path "versions.yml"                   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F'\t' 'NR==FNR {names[\$1]; next} \$1 in names' ${reps} ${clusters} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
