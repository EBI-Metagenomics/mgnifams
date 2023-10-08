process EXPORT_REPS {
    publishDir "${params.outDir}", mode: "copy"

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
    publishDir "${params.outDir}families/", mode: "copy"
    label "venv"

    input:
    path clust_tsv
    path fasta
    val mgyp

    output:
    path "known/${mgyp}_family.fa", optional: true, emit: known_ch
    path "unknown/${mgyp}_family.fa", optional: true, emit: unknown_ch

    script:
    """
    mkdir -p known
    mkdir -p unknown
    python3 ${params.scriptDir}family_rep_into_fasta.py ${clust_tsv} ${fasta} ${mgyp}_family.fa ${mgyp}
    """
}

process EXPORT_REPS_FA {
    publishDir "${params.outDir}", mode: "copy"
    label "seqtk"

    input:
    path fasta_files

    output:
    path "unknown_reps.fa"

    script:
    """
    for fasta in ${fasta_files.join(' ')}; do
        seqtk subseq \$fasta <(echo "\$(head -1 \$fasta | cut -c 2-)") >> unknown_reps.fa
    done
    """
}