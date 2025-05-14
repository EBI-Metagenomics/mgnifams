process PRESENT_DISCARDED_FAMILIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(discarded)

    output:
    tuple val(meta), path("${prefix}_discarded_mqc.csv"), emit: mqc
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<EOF > ${prefix}_discarded_mqc.csv
    # id: "discarded_clusters"
    # section_name: "Discarded clusters"
    # description: "Some of the initial clusters above minimum thresholds might not have created a family. If so, they are presented in a table below."
    # format: "csv"
    # plot_type: "table"
    Cluster Representative,Discard Reason,Value
    EOF

    cat ${discarded} >> ${prefix}_discarded_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_discarded_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
