process MAFFT {
    publishDir 'data/output/mafft', mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta}_mafft.fa"

    script:
    """
    mafft --anysymbol --thread 8 ${fasta} > ${fasta}_mafft.fa
    """
}