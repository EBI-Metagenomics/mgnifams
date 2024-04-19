process CREATE_CLUSTERS_PKL {
    publishDir "${params.outdir}/input/", mode: "copy"
    
    conda "${moduleDir}/environment.yml"
    
    input:
    path(clusters_tsv)

    output:
    path("clusters_bookkeeping_df.pkl"), emit: pkl
    path("log.txt")                    , emit: log

    script:
    """
    python3 ${params.scriptDir}/family/create_clusters_bookkeeping_df.py ${clusters_tsv}
    """
}

process REFINE_FAMILIES {
    publishDir "${params.outdir}/families/", mode: "copy"
    
    conda "${moduleDir}/environment.yml"
    
    input:
    path(anti_defence_proteins)
    path(clusters_pkl)
    path(families_tsv)
    path(fasta)
    path(discarded_clusters)
    path(converged_families)
    val(iteration)

    output:
    path("updated_refined_families.tsv")  , emit: tsv
    path("updated_mgnifams_dict.fa")      , emit: fa
    path("updated_discarded_clusters.txt"), emit: discarded
    path("seed_msa_sto/*")                , emit: seed_msa_sto
    path("msa_sto/*")                     , emit: msa_sto
    path("hmm/*")                         , emit: hmm
    path("domtblout/*")                   , emit: domtblout
    path("rf/*")                          , emit: rf
    path("updated_converged_families.txt"), emit: converged
    path("name_mapping.csv")              , emit: csv
    path("log.txt")                       , emit: log

    script:
    """
    python3 ${params.scriptDir}/family/refine_families.py ${anti_defence_proteins} ${clusters_pkl} ${families_tsv} ${fasta} ${discarded_clusters} ${converged_families} ${iteration}
    """
}

process EXTRACT_FIRST_STOCKHOLM_SEQUENCES {
    label "venv"

    input:
    path msa
    path ids
    val mode

    output:
    path "family_reps.fasta"

    script:
    """
    python3 ${params.scriptDir}/extract_first_stockholm_sequences.py ${msa} ${ids} ${mode} family_reps.fasta
    """
}
