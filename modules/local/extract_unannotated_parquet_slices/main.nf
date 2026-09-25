process EXTRACT_UNANNOTATED_PARQUET_SLICES {
    tag "${meta.id}_${chunk_index}"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861/data' :
        'community.wave.seqera.io/library/python_pyarrow:a33176f6cf91c593' }"

    input:
    tuple val(meta), path(sequences), path(pfam), val(chunk_index)
    val n_chunks
    val min_sequence_length

    output:
    tuple val(meta), path("${prefix}.faa"), emit: fa
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "chunk_${String.format('%06d', chunk_index)}"
    """
    extract_unannotated_parquet_slices.py \\
        --sequences "${sequences}" \\
        --pfam "${pfam}" \\
        --chunk_index ${chunk_index} \\
        --n_chunks ${n_chunks} \\
        --min_sequence_length ${min_sequence_length} \\
        --output_file "${prefix}.faa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "chunk_${String.format('%06d', chunk_index)}"
    """
    touch "${prefix}.faa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """
}
