process EXTRACT_CUDA_FAILED {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(scores, stageAs: "scores/*")

    output:
    tuple val(meta), path("cuda_failed_reps.fasta"), emit: fasta, optional: true
    path "versions.yml"                            , emit: versions

    script:
    """
    extract_cuda_failed.py \\
        --input_fasta ${fasta} \\
        --input_scores_folder scores \\
        --output_fasta cuda_failed_reps.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch cuda_failed_reps.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
