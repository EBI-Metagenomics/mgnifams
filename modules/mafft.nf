process MAFFT {
    publishDir "${params.outDir}mafft", mode: "copy", saveAs: { filename ->
        def newFilename = filename.replaceAll("family.fa_mafft", "msa")
        "${newFilename}"
    }
    label "mafft"

    input:
    path fasta

    output:
    path "${fasta}_mafft.fa"

    script:
    """
    mafft --anysymbol --thread 8 ${fasta} > ${fasta}_mafft.fa 2> /dev/null
    """
}