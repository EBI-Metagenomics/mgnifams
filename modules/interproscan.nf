process INTERPROSCAN {
    publishDir "${params.outDir}annotations/interproscan", mode: "copy"
    label "ips"

    input:
    path fasta

    output:
    path "${fasta}.tsv"

    script:
    """
    ${params.interproscan_dir}interproscan.sh \
    -cpu ${task.cpus} \
    -dp \
    -i ${fasta} \
    -f tsv
    """
    // --goterms \
    // -pa \
}