process EXPORT_MODEL_PFAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta ), path(hhr)

    output:
    tuple val(meta), path("mgnifam_pfams_${prefix}.csv"), emit: csv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = hhr.getExtension() == "gz" ? true : false
    def hhr_name = is_compressed ? hhr.getBaseName() : hhr
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${hhr} > ${hhr_name}
    fi

    awk '
        BEGIN { FS=";" }
        \$1 ~ /^>/ {
            gsub(/^>/, "", \$1)
            pfam = \$1
            short_desc = \$2
            long_desc = \$3
            gsub(/^[[:space:]]+/, "", short_desc)
            gsub(/^[[:space:]]+/, "", long_desc)
            if (!(pfam in seen)) {
                seen[pfam] = 1
                print pfam "\t" short_desc "\t" long_desc
            }
        }
    ' "${hhr_name}" > descriptions.tsv

    echo 'Descriptions mapping produced'

    export_pfams.py \\
        --id ${prefix} \\
        --hhr_file ${hhr_name} \\
        --descriptions descriptions.tsv \\
        --outfile mgnifam_pfams_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch mgnifam_pfams_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
