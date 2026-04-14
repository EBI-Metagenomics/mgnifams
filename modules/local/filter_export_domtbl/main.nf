process FILTER_EXPORT_DOMTBL {
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
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = domtbl.getExtension() == "gz" ? true : false
    def domtbl_name = is_compressed ? domtbl.getBaseName() : domtbl
    def header = (prefix == "pfam") ?
        "id,${prefix},name,e_value,score,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,acc" :
        "id,${prefix},e_value,score,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,acc"
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${domtbl} > ${domtbl_name}
    fi

    (
        echo "${header}"
        grep -v '^#' ${domtbl_name} | \\
            awk -v thr=${query_hmm_length_threshold} -v pfam=${prefix == "pfam" ? 1 : 0} '{
                qlen=\$6; env_from=\$20; env_to=\$21;
                env_len=env_to - env_from + 1;
                if (env_len >= thr * qlen) {
                    if (pfam) {
                        split(\$5, a, ".");   # split accession by "."
                        acc=a[1];             # keep only first part
                        printf "%s %s %s %s %s %s %s %s %s %s %s %s\\n",\$1,acc,\$4,\$7,\$8,\$16,\$17,\$18,\$19,\$20,\$21,\$22;
                    } else {
                        printf "%s %s %s %s %s %s %s %s %s %s %s\\n",\$1,\$4,\$7,\$8,\$16,\$17,\$18,\$19,\$20,\$21,\$22;
                    }
                }
            }' | tr ' ' ','
    ) > mgnifam_${prefix}s.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def header = (prefix == "pfam") ?
        "id,${prefix},name,e_value,score,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,acc" :
        "id,${prefix},e_value,score,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,acc"
    """
    echo "${header}" > mgnifam_${prefix}s.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
