process COUNT_SEED_MSA_SIZES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(seed_msa_sto, stageAs: 'seed_msa_sto/*')

    output:
    tuple val(meta), path("seed_sizes.csv"), emit: csv
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Seed size = unique sequence names per Stockholm file (skips #=G* markup, // and blank lines; multi-block safe)
    """
    set -euo pipefail

    echo "id,seed_size" > seed_sizes.csv
    find -L seed_msa_sto -type f -print0 | sort -z | xargs -0 -r awk '
        FNR == 1 { if (id != "") print id "," n; id = FILENAME; sub(/.*\\//, "", id); sub(/\\..*/, "", id); n = 0; split("", seen) }
        !/^#/ && !/^\\/\\// && NF && !seen[\$1]++ { n++ }
        END { if (id != "") print id "," n }
    ' >> seed_sizes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1 | cut -d, -f1)
    END_VERSIONS
    """

    stub:
    """
    echo "id,seed_size" > seed_sizes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1 | cut -d, -f1)
    END_VERSIONS
    """
}
