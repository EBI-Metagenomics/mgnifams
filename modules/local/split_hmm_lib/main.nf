process SPLIT_HMM_LIB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(lib)
    val chunk_size

    output:
    tuple val(meta), path("chunk_*.hmm.lib.gz"), emit: hmms
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Every NAME is a family id, i.e. the integer DB key: a positive int64, unique in the library
    """
    set -euo pipefail

    gzip -cdf "${lib}" | awk -v size=${chunk_size} '
        function fail(msg) { print "SPLIT_HMM_LIB: " msg > "/dev/stderr"; failed = 1; exit 1 }
        !in_record && NF == 0 { next }
        !in_record { in_record = 1; name = ""; out = sprintf("chunk_%06d.hmm.lib", int(records / size)) }
        \$1 == "NAME" {
            name = \$2 ""
            if (name !~ /^[1-9][0-9]*\$/ || length(name) > 19 || (length(name) == 19 && name > "9223372036854775807"))
                fail("NAME " name " is not a positive int64 family id")
            if (name in seen) fail("duplicate NAME " name)
            seen[name] = 1
        }
        { print > out }
        \$0 ~ /^\\/\\// {
            if (name == "") fail("HMM record " records + 1 " has no NAME")
            in_record = 0; records++
            if (records % size == 0) close(out)
        }
        END {
            if (failed) exit 1
            if (in_record) fail("truncated HMM record after " records " records")
            if (records == 0) fail("no HMM records in ${lib}")
        }
    '
    gzip chunk_*.hmm.lib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1 | cut -d, -f1)
        gzip: \$(gzip --version 2>&1 | head -n 1 | sed 's/^gzip //')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > chunk_000000.hmm.lib.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -n 1 | cut -d, -f1)
        gzip: \$(gzip --version 2>&1 | head -n 1 | sed 's/^gzip //')
    END_VERSIONS
    """
}
