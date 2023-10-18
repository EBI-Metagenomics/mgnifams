process HMMBUILD {
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll("_msa", "")
        "${newFilename}"
    }
    label "hmmer"

    input:
    path msa_folder

    output:
    path("hmm/*")

    script:
    """
    mkdir -p hmm
    for fasta in ${msa_folder}; do
        id=\$(basename \$fasta .fa)
        hmmbuild hmm/\${id}.hmm \$fasta
    done
    """
}

process HMMSCAN {
    publishDir "${params.outdir}hmm/scan", mode: "copy",
        pattern: "*_scan.tblout", saveAs: { filename ->
            def newFilename = filename.replaceAll(".fa_mafft.fa.hmm_scan", "_domains")
            "${newFilename}"
        }
    label "hmmer"

    input:
    path hmm
    path uniprot_fasta

    output:
    path "${hmm}.h3f"
    path "${hmm}.h3i"
    path "${hmm}.h3m"
    path "${hmm}.h3p"
    path "${hmm}_scan.tblout", emit: tblout_ch

    script:
    """
    hmmpress ${hmm} > /dev/null 2>&1
    hmmscan --domtblout ${hmm}_scan.tblout ${hmm} ${uniprot_fasta} > /dev/null
    """
}