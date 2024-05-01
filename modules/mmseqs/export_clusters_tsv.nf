process EXPORT_CLUSTERS_TSV {
    publishDir "${params.outDir}/mmseqs", mode: "copy"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    tuple val(meta1), path(mmseqs_DB)
    tuple val(meta2), path(mmseqs_clu)

    output:
    tuple val(meta2), path("${mmseqs_clu}.tsv"), emit: tsv

    script:
    """
    mmseqs \\
        createtsv \\
        ${mmseqs_DB}/${mmseqs_DB} \\
        ${mmseqs_DB}/${mmseqs_DB} \\
        ${mmseqs_clu}/${mmseqs_clu} \\
        ${mmseqs_clu}.tsv \\
        --threads ${task.cpus} \\
    """
}
