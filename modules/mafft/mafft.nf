process MAFFT {
    publishDir 'data/output', mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta}_mafft.fa"

    script:
    """
    mafft --thread 8 ${fasta} > ${fasta}_mafft.fa
    """
}