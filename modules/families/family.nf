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
    publishDir 'data/output', mode: 'copy'
    
    input:
    path clust_tsv
    path fasta
    val mgyp

    output:
    path "${clust_tsv}.family.fa"

    script:
    """
    python3 family_rep_into_fasta.py ${clust_tsv} ${fasta} ${clust_tsv}.family.fa ${mgyp}
    """
}