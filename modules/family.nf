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
    publishDir 'data/output/families', mode: 'copy'
    
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

process EXPORT_REPS_FA {
    publishDir 'data/output', mode: 'copy'

    input:
    path fasta_files

    output:
    path "reps.fa"

    script:
    """
    for fasta in ${fasta_files.join(' ')}; do
        seqtk subseq \$fasta <(echo "\$(head -1 \$fasta | cut -c 2-)") >> reps.fa
    done
    """
}