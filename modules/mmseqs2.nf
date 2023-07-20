process CREATEDB {
    publishDir 'data/output/mmseqs', mode: 'copy'

    input:
    path fasta

    output:
    path "DB/*"

    script:
    fastaName = "${fasta.getName().tokenize('.')[0]}"

    """
    mkdir DB
    mmseqs createdb ${fasta} DB/${fastaName}_DB > /dev/null 2>&1
    """
}

process LINCLUST {
    publishDir 'data/output/mmseqs', mode: 'copy'

    input:
    path mmseqs_DB

    output:
    path "clu/*"

    script:
    dbName = "${mmseqs_DB[0]}"

    """
    mkdir clu
    mmseqs linclust \
    ${dbName} \
    clu/${dbName}_clu \
    tmp \
    --min-seq-id ${params.seq_identity} \
    --cov-mode 1 \
    -c ${params.coverage} > /dev/null 2>&1
    """
}

process CREATETSV {
    publishDir 'data/output/mmseqs', mode: 'copy'
    
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
    ${cluName}.tsv > /dev/null 2>&1
    """
}

process CONVERT2FASTA {
    publishDir 'data/output/mmseqs', mode: 'copy'

    input:
    path mmseqs_DB
    path mmseqs_clu

    output:
    path "${cluName}_rep.fasta"

    script:
    dbName = "${mmseqs_DB[0]}"
    cluName = "${mmseqs_clu.baseName[0]}"

    """
    mmseqs createsubdb ${cluName} ${dbName} ${cluName}_rep
    mmseqs createsubdb ${cluName} ${dbName}_h ${cluName}_rep_h
    mmseqs convert2fasta ${cluName}_rep ${cluName}_rep.fasta
    """
}