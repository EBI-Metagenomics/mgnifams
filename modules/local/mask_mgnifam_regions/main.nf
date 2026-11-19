// MGnifam regions cut out of a protein FASTA, by hash buckets (bin/mask_mgnifam_regions.py).
// SPLIT_INTO_BUCKETS writes every bucket file (empty or not), so MASK_BUCKET gets one FASTA and one region file per bucket.

process SPLIT_INTO_BUCKETS {
    tag "${meta.id}_${kind}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data':
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(inputs, stageAs: 'input/*') // FASTA files, or refined_families-style region TSVs
    val kind                                          // 'fasta' or 'regions'
    val buckets

    output:
    tuple val(meta), path("buckets/*"), emit: buckets // buckets/000.fa ... or buckets/000.tsv ...
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mask_mgnifam_regions.py split --kind ${kind} --buckets ${buckets} --out buckets ${inputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def ext = kind == 'fasta' ? 'fa' : 'tsv'
    """
    mkdir buckets
    for i in \$(seq 0 \$(( ${buckets} - 1 ))); do touch buckets/\$(printf '%03d' \$i).${ext}; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}

process MASK_BUCKET {
    tag "${meta.id}_${bucket}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data':
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), val(bucket), path(fasta), path(regions)
    val min_sequence_length

    output:
    tuple val(meta), path("${bucket}.masked.fa"), emit: fa
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mask_mgnifam_regions.py mask \\
        --fasta ${fasta} \\
        --regions ${regions} \\
        --min_sequence_length ${min_sequence_length} \\
        --out ${bucket}.masked.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${bucket}.masked.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
