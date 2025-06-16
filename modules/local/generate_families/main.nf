process GENERATE_FAMILIES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e26ce00244b6cc71edc451bb368f0e4e6ad2cf497f152337177ea0bff7d5cba1/data' :
        'community.wave.seqera.io/library/python_pip_biopython_numpy_pruned:dd4bec78b7f08a4e' }"
    
    input:
    tuple val(meta) , path(clusters_chunk)
    tuple val(meta2), path(mgnifams_fasta)
    val(discard_min_rep_length)
    val(discard_max_rep_length)
    val(discard_min_starting_membership)
    val(max_seq_identity)
    val(max_seed_seqs)
    val(max_gap_occupancy)
    val(recruit_evalue_cutoff)
    val(recruit_hit_length_percentage)

    output:
    tuple val(meta), path("seed_msa_sto/*")       , emit: seed_msa_sto
    tuple val(meta), path("full_msa_sto/*")       , emit: full_msa_sto
    tuple val(meta), path("hmm/*")                , emit: hmm
    tuple val(meta), path("rf/*")                 , emit: rf
    tuple val(meta), path("refined_families/*")   , emit: tsv
    tuple val(meta), path("discarded_clusters/*") , emit: discarded
    tuple val(meta), path("successful_clusters/*"), emit: successful
    tuple val(meta), path("converged_families/*") , emit: converged
    tuple val(meta), path("family_metadata/*")    , emit: metadata
    tuple val(meta), path("logs/*")               , emit: logs
    tuple val(meta), path("family_reps/*")        , emit: fasta
    path("versions.yml")                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_families.py \\
        --clusters_chunk ${clusters_chunk} \\
        --fasta_file ${mgnifams_fasta} \\
        --cpus ${task.cpus} \\
        --chunk_num ${prefix} \\
        --discard_min_rep_length ${discard_min_rep_length} \\
        --discard_max_rep_length ${discard_max_rep_length} \\
        --discard_min_starting_membership ${discard_min_starting_membership} \\
        --max_seq_identity ${max_seq_identity} \\
        --max_seed_seqs ${max_seed_seqs} \\
        --max_gap_occupancy ${max_gap_occupancy} \\
        --recruit_evalue_cutoff ${recruit_evalue_cutoff} \\
        --recruit_hit_length_percentage ${recruit_hit_length_percentage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
        pyfamsa: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfamsa'))")
        pyhmmer: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyhmmer'))")
        pytrimal: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pytrimal'))")
        numpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('numpy'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    mkdir seed_msa_sto
    cp ${prefix}.txt seed_msa_sto/
    mkdir full_msa_sto
    cp ${prefix}.txt full_msa_sto/
    mkdir hmm
    cp ${prefix}.txt hmm/
    mkdir rf
    cp ${prefix}.txt rf/
    mkdir domtblout
    cp ${prefix}.txt domtblout/
    mkdir refined_families
    cp ${prefix}.txt refined_families/
    mkdir discarded_clusters
    cp ${prefix}.txt discarded_clusters/
    mkdir successful_clusters
    cp ${prefix}.txt successful_clusters/
    mkdir converged_families
    cp ${prefix}.txt converged_families/
    mkdir family_metadata
    cp ${prefix}.txt family_metadata/
    mkdir logs
    cp ${prefix}.txt logs/
    mkdir family_reps
    mv ${prefix}.txt family_reps/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
        pyfastx: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfastx'))")
        pyfamsa: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfamsa'))")
        pyhmmer: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyhmmer'))")
        pytrimal: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pytrimal'))")
        numpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('numpy'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
