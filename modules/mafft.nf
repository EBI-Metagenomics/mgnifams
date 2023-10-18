process MAFFT {
    publishDir "${params.outdir}", mode: 'copy'
    label "mafft"

    input:
    path non_singletons_folder

    output:
    path("mafft/*")

    script:
    """
    mkdir -p mafft
    for fasta in ${non_singletons_folder}; do
        id=\$(basename \$fasta .fasta)
        mafft --anysymbol --thread 8 \$fasta > mafft/\${id}_msa.fa 2> /dev/null
    done
    """
}
