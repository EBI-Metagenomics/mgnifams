process CREATE_DB {
    input:
    path fasta

    output:
    path "mmseqs_db/*"

    script:
    """
    mkdir mmseqs_db
    mmseqs createdb ${fasta} mmseqs_db/${fasta}_db > /dev/null 2>&1
    """
}

process LINCLUST {
    input:
    path mmseqsDB

    output:
    path "mmseqs_clus/*"

    script:
    """
    mkdir mmseqs_clus
    mmseqs linclust ${mmseqsDB.baseName[1]} mmseqs_clus/${mmseqsDB.baseName[0]}_clu mmseqs_clus/tmp --min-seq-id 0.3 --cov-mode 1 -c 0.8 > /dev/null 2>&1
    """
}

process CREATE_TSV {
    publishDir 'data/output', mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll(".fa", "")
        "${newFilename}"
    }
    
    input:
    path mmseqsDB
    path mmseqsCLU

    output:
    path "${mmseqsCLU.baseName[0]}.tsv"

    script:
    """
    mmseqs createtsv ${mmseqsDB.baseName[1]} ${mmseqsDB.baseName[1]} ${mmseqsCLU.baseName[0]} ${mmseqsCLU.baseName[0]}.tsv > /dev/null 2>&1
    """
}