process PARSE_FAMILIES {
    publishDir "${params.outDir}families/", mode: "copy"
    label "venv"

    input:
    path fasta
    path clust_tsv

    output:
    path "rep_names.txt"   , emit: reps_ids
    path "reps.fa"         , emit: reps_fasta
    path "singletons.txt"  , optional: true, emit: singleton_ids
    path "non_singletons/*", optional: true, emit: non_singletons_ch

    script:
    """
    python3 ${params.scriptDir}parse_families.py ${fasta} ${clust_tsv}
    """
}
