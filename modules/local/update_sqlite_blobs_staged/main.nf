process UPDATE_SQLITE_BLOBS_STAGED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/313e1c18a344323886cf97a151ab66d81c1a146fb129558cb9382b69a72d5532/data' :
        'community.wave.seqera.io/library/python:b1b4b1f458c605bb' }"

    // Blob inputs are staged per column (each may be []); files are <family_id>.<ext>, found recursively
    input:
    tuple val(meta), path(db, stageAs: 'input/db.sqlite3')
    path cifs   , stageAs: 'cif/*'
    path s4pred , stageAs: 's4pred/*'
    path tm     , stageAs: 'tm/*'
    path biomes , stageAs: 'biome/*'
    path domains, stageAs: 'domain/*'
    path ids                     // expected family ids, one per line
    val info                     // Map for the update_info table; tm_computed/biome_computed make those blobs required

    output:
    tuple val(meta), path("${prefix}_update.sqlite3"), emit: db
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def required = ['cif', 's4pred', 'domain'] + (info.tm_computed ? ['tm'] : []) + (info.biome_computed ? ['biome'] : [])
    def info_args = info.collect { k, v -> "'${k}=${v}'" }.join(' ')
    """
    mkdir -p cif s4pred tm biome domain
    cp -L "${db}" "${prefix}_update.sqlite3"

    update_sqlite_blobs_staged.py \\
        --db "${prefix}_update.sqlite3" \\
        --cif_dir cif \\
        --s4pred_dir s4pred \\
        --tm_dir tm \\
        --biome_dir biome \\
        --domain_dir domain \\
        --ids "${ids}" \\
        --required ${required.join(',')} \\
        --info ${info_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -L "${db}" "${prefix}_update.sqlite3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
