process CREATEDB {
    publishDir "${params.outDir}mmseqs", mode: "copy"
    label "mmseqs"

    input:
    path fasta

    output:
    path "${fastaName}_DB/*"

    script:
    fastaName = "${fasta.getName().tokenize('.')[0]}"

    """
    mkdir ${fastaName}_DB
    mmseqs createdb ${fasta} ${fastaName}_DB/${fastaName}_DB > /dev/null 2>&1
    """
}

process LINCLUST {
    publishDir "${params.outDir}mmseqs", mode: "copy"
    label "mmseqs"

    input:
    path mmseqs_DB

    output:
    path "${dbName}_clu/*"

    script:
    dbName = "${mmseqs_DB[0]}"

    """
    mkdir ${dbName}_clu
    mmseqs linclust \
    ${dbName} \
    ${dbName}_clu/${dbName}_clu \
    tmp \
    --min-seq-id ${params.seq_identity} \
    --cov-mode 1 \
    -c ${params.coverage} \
    --threads ${task.cpus} > /dev/null 2>&1
    """
}

process CREATETSV {
    publishDir "${params.outDir}mmseqs", mode: "copy"
    label "mmseqs"

    input:
    path mmseqs_DB
    path mmseqs_clu

    output:
    path "${cluName}.tsv"

    script:
    dbName = "${mmseqs_DB[0]}"
    cluName = "${mmseqs_clu.baseName[0]}"

    """
    mmseqs createtsv \
    ${dbName} \
    ${dbName} \
    ${cluName} \
    ${cluName}.tsv \
    --threads ${task.cpus} > /dev/null 2>&1
    """
}

process CONVERT2FASTA {
    publishDir "${params.outDir}mmseqs", mode: "copy"
    label "mmseqs"

    input:
    path mmseqs_DB
    path mmseqs_clu

    output:
    path "${cluName}_rep.fasta"

    script:
    dbName = "${mmseqs_DB[0]}"
    cluName = "${mmseqs_clu.baseName[0]}"

    """
    mmseqs createsubdb ${cluName} ${dbName} ${cluName}_rep > /dev/null 2>&1
    mmseqs createsubdb ${cluName} ${dbName}_h ${cluName}_rep_h > /dev/null 2>&1
    mmseqs convert2fasta ${cluName}_rep ${cluName}_rep.fasta > /dev/null 2>&1
    """
}

process CLUSTERUPDATE {
    publishDir "${params.outDir}mmseqs", mode: "copy"
    label "mmseqs"
    
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