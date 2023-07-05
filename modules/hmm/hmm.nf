process HMMBUILD {
    publishDir 'data/output/hmm', mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta}.hmm"

    script:
    """
    hmmbuild ${fasta}.hmm ${fasta}
    """
}

process HMMSCAN {
    publishDir 'data/output/hmm', pattern: '*_scan.txt', mode: 'copy'

    input:
    path hmm
    path uniprot_fasta

    output:
    path "${hmm}.h3f"
    path "${hmm}.h3i"
    path "${hmm}.h3m"
    path "${hmm}.h3p"
    path "${hmm}_scan.txt"

    script:
    """
    hmmpress ${hmm}
    hmmscan ${hmm} ${uniprot_fasta} > ${hmm}_scan.txt
    """
}