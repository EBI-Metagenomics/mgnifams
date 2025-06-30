process EXTRACT_ESMFOLD_SCORES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(scores)

    output:
    tuple val(meta), path("${meta.id}_scores.csv"), emit: csv
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    csv_file="${meta.id}_scores.csv"
    echo "family_id,rep_length,plddt,ptm" > "\$csv_file"
    grep -E 'pLDDT' *_scores.txt | awk -F' |, ' '{gsub("_", "/", \$11); print \$11 "," \$14 "," \$16 "," \$18}' >> "\$csv_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clustering_distribution_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
