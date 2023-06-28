process EXPORT_REPS {
    // TODO, for non 1-member clusters
}

process CREATE_FAMILY_FA {
    publishDir 'data/output', mode: 'copy'
    
    input:
    path clusters
    path fasta
    val mgyp

    output:
    path "${clusters}.family.fa"

    script:
    """
    python3 family_rep_into_fasta.py ${clusters} ${fasta} ${clusters}.family.fa ${mgyp}
    """
}