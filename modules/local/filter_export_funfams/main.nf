process FILTER_EXPORT_FUNFAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta ), path(domtbl)
    val(query_hmm_length_threshold)

    output:
    tuple val(meta), path("mgnifam_funfams.csv"), emit: csv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = domtbl.getExtension() == "gz" ? true : false
    def domtbl_name = is_compressed ? domtbl.getBaseName() : domtbl
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${domtbl} > ${domtbl_name}
    fi

    (
        echo "id,funfam,e_value,score,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,acc"
        grep -v '^#' ${domtbl_name} | \
            awk '{
                \$1=\$1; gsub(/ +/, " ");
                qlen=\$6; env_from=\$20; env_to=\$21;
                env_len=env_to - env_from + 1;
                if (env_len >= ${query_hmm_length_threshold} * qlen) print
            }' | \
            cut -d' ' -f1,4,7,8,16,17,18,19,20,21,22 | \
            tr ' ' ','
    ) > mgnifam_funfams.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    """
    touch mgnifam_funfams.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
