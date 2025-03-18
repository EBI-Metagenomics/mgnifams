process GENERATE_FAMILIES {
    tag "$meta.id"
    label 'process_high'
    maxForks 30 // TODO from config

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/microbiome-informatics/mgnifams:3.2.0' :
        'quay.io/microbiome-informatics/mgnifams:3.2.0' }"
    
    input:
    tuple val(meta) , path(clusters_chunk)
    tuple val(meta2), path(mgnifams_fasta)

    output:
    tuple val(meta), path("seed_msa_sto/*")       , emit: seed_msa_sto
    tuple val(meta), path("msa_sto/*")            , emit: msa_sto
    tuple val(meta), path("hmm/*")                , emit: hmm
    tuple val(meta), path("rf/*")                 , emit: rf
    tuple val(meta), path("domtblout/*")          , emit: domtblout
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
        --chunk_num ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
        pyfastx: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfastx'))")
        pyfamsa: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfamsa'))")
        pyhmmer: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyhmmer'))")
        pytrimal: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pytrimal'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    mkdir seed_msa_sto
    cp ${prefix}.txt seed_msa_sto/
    mkdir msa_sto
    cp ${prefix}.txt msa_sto/
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
    mv ${prefix}.txt logs/
    mkdir family_reps
    mv ${prefix}.fasta family_reps/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
        pyfastx: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfastx'))")
        pyfamsa: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfamsa'))")
        pyhmmer: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyhmmer'))")
        pytrimal: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pytrimal'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
