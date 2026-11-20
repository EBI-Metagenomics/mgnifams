process FORMAT_CLANS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    input:
    tuple val(meta), path(clans), path(clan_edges)

    output:
    tuple val(meta), path("clans.tsv.gz")     , emit: clans
    tuple val(meta), path("clan_edges.tsv.gz"), emit: clan_edges
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | sed '1!d; s/gzip //'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    // The FTP clan files: compute_clans CSVs as tab-separated gzip (no field holds a comma); -n keeps md5s stable
    """
    tr , '\\t' < ${clans} | gzip -n > clans.tsv.gz
    tr , '\\t' < ${clan_edges} | gzip -n > clan_edges.tsv.gz
    """

    stub:
    """
    echo | gzip -n > clans.tsv.gz
    echo | gzip -n > clan_edges.tsv.gz
    """
}
