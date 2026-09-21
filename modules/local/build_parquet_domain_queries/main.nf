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

    # Pfam accession (no version) -> description, as in the sequence_explorer_pfam mapping.
    # Python, not gzip: the python_pyarrow container has no gzip binary.
    python - "${pfam_hmms}" > pfam_mapping.tsv <<-'END_MAPPING'
    import gzip, sys

    path = sys.argv[1]
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\\x1f\\x8b'
    acc = ''
    with (gzip.open(path, 'rt') if gzipped else open(path)) as fh:
        for line in fh:
            if line.startswith('ACC'):
                acc = line.split()[1].split('.')[0]
            elif line.startswith('DESC'):
                print(acc + '\\t' + line.split(None, 1)[1].rstrip('\\n'))
    END_MAPPING

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
