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

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # MMseqs2 TSV row order changes between runs; sort by rep, then member with the rep first,
    # so mgnifam numbers families (and FAMSA reads members) in the same order every run
    awk -F'\t' -v OFS='\t' 'NR==FNR {names[\$1]; next} \$1 in names {print \$1, (\$1 != \$2), \$2}' ${reps} ${clusters} \\
        | LC_ALL=C sort -t "\$(printf '\\t')" -k1,1 -k2,2 -k3,3 \\
        | cut -f1,3 > ${prefix}.tsv

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
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
