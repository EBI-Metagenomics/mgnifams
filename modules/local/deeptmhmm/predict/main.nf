process DEEPTMHMM_PREDICT {
    tag "$meta.id"
    label 'process_single'
    label 'process_gpu'
    stageInMode 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://interpro/deeptmhmm:1.0' :
        'docker.io/interpro/deeptmhmm:1.0' }"

    input:
    tuple val(meta), path(fasta)
    path tmhmm_dir

    output:
    tuple val(meta), path("biolib_results")
    tuple val(meta), path("biolib_results/TMRs.gff3")                 , emit: gff3
    tuple val(meta), path("biolib_results/predicted_topologies.3line"), emit: line3
    tuple val(meta), path("biolib_results/deeptmhmm_results.md")      , emit: md
    tuple val(meta), path("biolib_results/*_probs.csv")               , optional: true, emit: csv
    tuple val(meta), path("biolib_results/plot.png")                  , optional: true, emit: png
    path "versions.yml"                                               , emit: versions

    script:
    """
    # deeptmhmm has a hard coded assumption it is being run within its dir
    cd ${tmhmm_dir}
    python3 predict.py \\
        --fasta ../${fasta} \\
        --output-dir ../biolib_results
    cd ..
    rm -r biolib_results/embeddings
    chmod -R 777 biolib_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        numpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('numpy'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
        torch: \$(python -c "import importlib.metadata; print(importlib.metadata.version('torch'))")
        tqdm: \$(python -c "import importlib.metadata; print(importlib.metadata.version('tqdm'))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir biolib_results
    touch biolib_results/TMRs.gff3
    touch biolib_results/predicted_topologies.3line
    touch biolib_results/deeptmhmm_results.md
    touch biolib_results/MX_probs.csv
    touch biolib_results/plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        numpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('numpy'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
        torch: \$(python -c "import importlib.metadata; print(importlib.metadata.version('torch'))")
        tqdm: \$(python -c "import importlib.metadata; print(importlib.metadata.version('tqdm'))")
    END_VERSIONS
    """
}