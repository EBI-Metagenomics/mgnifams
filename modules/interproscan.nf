process INTERPROSCAN {
    publishDir "${params.outDir}annotations/interproscan", mode: "copy"

    input:
    path fasta

    output:
    path "${fasta}.tsv"

    script:
    """
    ${params.interproscan_dir}interproscan.sh -i ${fasta} -f tsv -dp
    """
}