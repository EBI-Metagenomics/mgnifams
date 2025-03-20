process PRESENT_FAMILY_METADATA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(metadata)

    output:
    tuple val(meta), path("${prefix}_metadata_mqc.csv"), emit: mqc
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<EOF > ${prefix}_metadata_mqc.csv
    # id: "family_metadata"
    # section_name: "Family metadata"
    # description: "Generated family metadata."
    # format: "csv"
    # plot_type: "table"
    Family Id,Size,Representative Id,Region,Representative Length,Sequence
    EOF

    cat ${metadata} >> ${prefix}_metadata_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metadata_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
