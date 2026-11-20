process CLASSIFY_MULTI_HITS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(full_msas, stageAs: "full_msa/*")

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // The trailing "# total ..." line is per chunk: dropped, the header is kept once when the chunks are collected
    """
    set -o pipefail
    classify_multi_hits.py full_msa/* | grep -v '^#' > ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'family\\tmodel_length\\tregions\\tmgyps\\tfull_length_regions\\tsingle_full\\tsingle_partial\\tmulti_hit_mgyps\\tsplit\\tsplit_full\\ttandem\\tmixed\\tpartial_repeat\\tcopies\\tfull_copies\\n' > ${prefix}.tsv
    """
}
