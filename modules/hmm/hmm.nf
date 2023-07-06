process HMMBUILD {
    publishDir 'data/output/hmm/build', mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll(".fa_mafft.fa", "")
        "${newFilename}"
    }
    
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
    publishDir 'data/output/hmm/scan', mode: 'copy',
        pattern: '*_scan.tblout', saveAs: { filename ->
            def newFilename = filename.replaceAll(".fa_mafft.fa.hmm_scan", "_domains")
            "${newFilename}"
        }

    input:
    path hmm
    path uniprot_fasta

    output:
    path "${hmm}.h3f"
    path "${hmm}.h3i"
    path "${hmm}.h3m"
    path "${hmm}.h3p"
    path "${hmm}_scan.tblout"

    script:
    """
    hmmpress ${hmm}
    hmmscan --domtblout ${hmm}_scan.tblout ${hmm} ${uniprot_fasta} > /dev/null
    """
}