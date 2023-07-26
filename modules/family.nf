process EXPORT_REPS {
    publishDir 'data/output', mode: 'copy'

    input:
    path clust_tsv

    output:
    path "rep_names.txt"

    script:
    """
    awk -F'\t' '{print \$1}' ${clust_tsv} | sort -u > rep_names.txt
    """
}

process CREATE_FAMILY_FA {
    publishDir 'data/output/families/all', mode: 'copy'
    
    input:
    path clust_tsv
    path fasta
    val mgyp

    output:
    path "${mgyp}_family.fa"

    script:
    """
    python3 ${baseDir}/bin/family_rep_into_fasta.py ${clust_tsv} ${fasta} ${mgyp}_family.fa ${mgyp}
    """
}

process KEEP_UNKNOWN {
    publishDir 'data/output/families/unknown', mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta.baseName}_u.fa", optional: true

    script:
    """
    python3 ${baseDir}/bin/keep_unknown_family.py ${fasta} ${fasta.baseName}_u.fa
    """
}