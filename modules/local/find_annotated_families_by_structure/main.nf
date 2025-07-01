process FIND_ANNOTATED_FAMILIES_BY_STRUCTURE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(m8s)

    output:
    tuple val(meta), path("annotated_families.txt"), emit: annotated_families, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    awk '{print \$1}' ${m8s} | cut -d'-' -f 1 | sort | uniq > annotated_families.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    """
    touch annotated_families.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
