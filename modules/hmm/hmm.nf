process HMMBUILD {
    publishDir 'data/output', mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta}.hmm"

    script:
    """
    hmmbuild ${fasta}.hmm ${fasta}
    """
}

// process HMMPRESS {}
// process HMMSCAN {}