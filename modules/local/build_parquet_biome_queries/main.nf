process BUILD_PARQUET_BIOME_QUERIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e415e86f361373d5bf62a587c79daada9fa797ce076c764f1ccadb32ba2c3861/data' :
        'community.wave.seqera.io/library/python_pyarrow:a33176f6cf91c593' }"

    input:
    tuple val(meta), path(families), path(clusters)

    output:
    tuple val(meta), path("query_results/*"), emit: res
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_parquet_biome_queries.py \\
        --families ${families} \\
        --clusters "${clusters}" \\
        --output_dir query_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir query_results
    cat ${families} | cut -f1 | sort -u | while read -r family; do touch "query_results/\${family}.tsv"; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """
}
