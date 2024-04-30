process CREATETSV {
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

process CLUSTERUPDATE {
    publishDir "${params.outDir}/mmseqs", mode: "copy"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    path mmseqs_DB_old
    path mmseqs_DB_new
    path mmseqs_DB_old_clu

    output:
    path "${new_dbName}_update/DB/*"
    path "${new_dbName}_update/clu/*"

    script:
    old_dbName = "${mmseqs_DB_old.baseName[0]}"
    new_dbName = "${mmseqs_DB_new.baseName[0]}"
    old_cluName = "${mmseqs_DB_old_clu.baseName[0]}"

    """
    mkdir ${new_dbName}_update
    mkdir ${new_dbName}_update/DB
    mkdir ${new_dbName}_update/clu
    mmseqs clusterupdate \
    ${old_dbName} \
    ${new_dbName} \
    ${old_cluName} \
    ${new_dbName}_update/DB/DB_new_updated \
    ${new_dbName}_update/clu/DB_update_clu \
    tmp \
    --threads ${task.cpus} > /dev/null 2>&1
    """
}
