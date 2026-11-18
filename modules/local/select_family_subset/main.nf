process SELECT_FAMILY_SUBSET {
    tag "$meta.id"
    label 'process_single'

    // ponytail: stdlib only; reuses the PARSE_CIF image instead of pulling another
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/aefa332c17483d2a371e10103c221676ad7deec3460ce9836d8bd0ed79c8d1d2/data' :
        'community.wave.seqera.io/library/python_pip_biopython_numpy:e51160118726626c' }"

    input:
    tuple val(meta), path(families, stageAs: 'families/*'), path(reps, stageAs: 'reps/*'), path(structures, stageAs: 'structures_in')

    output:
    tuple val(meta), path("family_ids.fasta"), emit: fasta
    tuple val(meta), path("pdb/*.pdb")       , emit: pdb, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def status = task.ext.status ?: 'not_updated'
    """
    select_family_subset.py \\
        --families ${families} \\
        --reps ${reps} \\
        --structures ${structures} \\
        --status ${status} \\
        --out_fasta family_ids.fasta \\
        --out_structures pdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    // The real script: stdlib only and fast, so stub runs select the same families as real ones
    stub:
    def status = task.ext.status ?: 'not_updated'
    """
    select_family_subset.py \\
        --families ${families} \\
        --reps ${reps} \\
        --structures ${structures} \\
        --status ${status} \\
        --out_fasta family_ids.fasta \\
        --out_structures pdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
