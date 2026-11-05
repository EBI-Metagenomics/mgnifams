process SPLIT_FAMILY_MEMBERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    // family<TAB>member TSV whose families are contiguous runs (refined_families.tsv of POOL_UPDATED_FAMILIES)
    input:
    tuple val(meta), path(tsv)
    val n_chunks

    output:
    tuple val(meta), path("chunks/*.tsv"), emit: chunks
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    mkdir chunks

    # Streams once, in constant memory (plus one entry per family id): about n_chunks chunks of equal bytes, a
    # chunk only ever starting at a family's first line, so a family never spans two chunks.
    # \$1 "" compares as strings: awk compares number-like fields numerically (9 == 09)
    awk -F '\\t' -v n=${n_chunks} -v size="\$(stat -L -c %s ${tsv})" '
        \$1 "" != prev {
            if (\$1 in seen) { print "${tsv}: family " \$1 " is not contiguous (line " NR ")" > "/dev/stderr"; exit 1 }
            seen[\$1] = 1
            prev = \$1 ""
            k = int(bytes * n / size)
            name = sprintf("chunks/families_%05d.tsv", k)
            if (name != out) { if (out) close(out); out = name }
        }
        { print > out; bytes += length(\$0) + 1 }
    ' ${tsv}
    ls chunks/*.tsv > /dev/null # an empty input fails here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir chunks
    cp -L ${tsv} chunks/families_00000.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
