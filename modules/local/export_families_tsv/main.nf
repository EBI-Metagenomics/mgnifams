process EXPORT_FAMILIES_TSV {
    tag "$meta.id"
    label 'process_single'

    // ponytail: stdlib only; reuses the PARSE_CIF image instead of pulling another
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/aefa332c17483d2a371e10103c221676ad7deec3460ce9836d8bd0ed79c8d1d2/data' :
        'community.wave.seqera.io/library/python_pip_biopython_numpy:e51160118726626c' }"

    input:
    tuple val(meta), path(previous, stageAs: 'previous/*'), path(delta), path(metadata)
    val mgnifams_release
    val mgnify_proteins_release

    output:
    tuple val(meta), path("families.tsv.gz"), emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export_families_tsv.py \\
        --previous ${previous} \\
        --delta ${delta} \\
        --metadata ${metadata} \\
        --release ${mgnifams_release} \\
        --proteins_release ${mgnify_proteins_release} \\
        --output families.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > families.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
