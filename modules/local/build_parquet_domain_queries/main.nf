process BUILD_PARQUET_DOMAIN_QUERIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyarrow:f6040e789f1cc7b8' :
        'community.wave.seqera.io/library/python_pyarrow:a33176f6cf91c593' }"

    input:
    tuple val(meta), path(refined_families), path(pfam)
    path pfam_hmms

    output:
    tuple val(meta), path("query_results/*")  , emit: res
    tuple val(meta), path("pfam_mapping.tsv") , emit: pfam_mapping
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    build_parquet_domain_queries.py \\
        --refined_families "${refined_families}" \\
        --pfam "${pfam}" \\
        --output_dir query_results

    # Pfam accession (no version) -> description, as in the sequence_explorer_pfam mapping
    gzip -cdf "${pfam_hmms}" | awk '/^ACC/ { split(\$2, a, "."); acc = a[1] } /^DESC/ { sub(/^DESC +/, ""); print acc "\\t" \$0 }' > pfam_mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir query_results
    cut -f1 "${refined_families}" | sort -u | while read -r family; do touch "query_results/\${family}.tsv"; done
    touch pfam_mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """
}
