process EXPORT_FOLDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta ), path(m8)

    output:
    tuple val(meta), path("mgnifam_folds.csv"), emit: csv
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = m8.getExtension() == "gz" ? true : false
    def m8_name = is_compressed ? m8.getBaseName() : m8
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${m8} > ${m8_name}
    fi

    (
        echo "id,fold,aligned_length,q_start,q_end,t_start,t_end,e_value"
        cut -f1,2,4,7,8,9,10,11 ${m8_name} | tr '\t' ','
    ) > mgnifam_folds.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    """
    touch mgnifam_folds.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
