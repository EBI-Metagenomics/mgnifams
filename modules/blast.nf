process MAKEBLASTDB {
    label "blast"

    input:
    path fasta

    output:
    path "${fasta}_DB/*"

    script:
    """
    makeblastdb \
    -in ${fasta} \
    -dbtype prot \
    -out ${fasta}_DB/${fasta}_DB
    """
}

process BLASTP {
    publishDir "${params.outDir}annotations/blast", mode: "copy"
    label "blast"

    input:
    path fasta
    path library

    output:
    path "${fasta}.csv"

    script:
    dbName = "${library.baseName[0]}"

    """
    blastp \
    -db ${dbName} \
    -query ${fasta} \
    -out ${fasta}.csv \
    -num_threads ${task.cpus} \
    -outfmt 10
    """
    // -outfmt 10: CSV
    // -evalue 1e-5 TODO how much?
}