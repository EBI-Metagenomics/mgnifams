process PARSE_CIF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/aefa332c17483d2a371e10103c221676ad7deec3460ce9836d8bd0ed79c8d1d2/data' :
        'community.wave.seqera.io/library/python_pip_biopython_numpy:e51160118726626c' }"

    input:
    tuple val(meta), path(pdb, stageAs: "pdb_folder/*")

    output:
    tuple val(meta), path("*.cif"), emit: cif
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """    
    for file in pdb_folder/*; do
        name=\$(basename \$file .pdb)
        parse_cif.py \\
            --pdb_file \$file \\
            --output_file \$name.cif
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        numpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('numpy'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clustering_distribution_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        numpy: \$(python -c "import importlib.metadata; print(importlib.metadata.version('numpy'))")
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
