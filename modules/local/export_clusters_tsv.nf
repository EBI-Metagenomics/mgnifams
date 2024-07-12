process EXPORT_CLUSTERS_TSV {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    tuple val(meta1), path(mmseqs_DB)
    tuple val(meta2), path(mmseqs_clu)

    output:
    tuple val(meta2), path("${mmseqs_clu}.tsv"), emit: tsv
    env num_sequences                          , emit: num_sequences
    path "versions.yml"                        , topic: 'versions'

    script:
    """
    mmseqs \\
        createtsv \\
        ${mmseqs_DB}/${mmseqs_DB} \\
        ${mmseqs_DB}/${mmseqs_DB} \\
        ${mmseqs_clu}/${mmseqs_clu} \\
        ${mmseqs_clu}.tsv \\
        --threads ${task.cpus}
    
    num_sequences=\$(wc -l < "${mmseqs_clu}.tsv")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
