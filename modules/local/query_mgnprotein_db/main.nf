process QUERY_MGNPROTEIN_DB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56121541d6b4f7082058be7bf67ed33c9e6d82dd57015ca431f632b696117179/data' :
        'community.wave.seqera.io/library/psycopg2:2.9.10--a5479766f9fb64b2' }"

    input:
    tuple val(meta), path(config_file), path(family_proteins_file)

    output:
    tuple val(meta), path("query_results/*")  , emit: res
    tuple val(meta), path("biome_mapping.tsv"), emit: biome_mapping
    tuple val(meta), path("pfam_mapping.tsv") , emit: pfam_mapping
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    query_mgnprotein_db.py \\
        --mgnprotein_db_config_file ${config_file} \\
        --family_proteins_file ${family_proteins_file} \\
        --output_dir query_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        psycopg2: \$(python -c "import importlib.metadata; print(importlib.metadata.version('psycopg2'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir query_results
    touch query_results/test.tsv
    touch biome_mapping.tsv
    touch pfam_mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        psycopg2: \$(python -c "import importlib.metadata; print(importlib.metadata.version('psycopg2'))")
    END_VERSIONS
    """
}
