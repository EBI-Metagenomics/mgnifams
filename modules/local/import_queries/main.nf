process IMPORT_QUERIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d41320ac1ca5f0a626982296d23ded50376a966b9e8240aa50dba2014a805bf5/data':
        'community.wave.seqera.io/library/sqlite:3.48.0--48957425ca78aa09' }"

    input:
    tuple val(meta) , path(pipeline_results)
    tuple val(meta2), path(db, stageAs: 'input/db.sqlite3')

    output:
    tuple val(meta), path("${prefix}.sqlite3"), emit: db
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def required = task.ext.args ?: '' // child CSVs that must be present, e.g. 'mgnifam_pfams.csv mgnifam_folds.csv'
    """
    set -euo pipefail

    for f in mgnifam.csv ${required}; do
        [ -f "\$f" ] || { echo "IMPORT_QUERIES: required \$f is missing" >&2; exit 1; }
    done

    # Work on a copy, so the INIT_SQLITE output stays pristine (-resume safe)
    cp -L "${db}" "${prefix}.sqlite3"

    cols() { sqlite3 "${prefix}.sqlite3" "SELECT group_concat(name, ',') FROM pragma_table_info('\$1') WHERE name != 'id'"; }
    nullif() { local IFS=,; local out=(); for c in \$1; do out+=("NULLIF(\$c, '')"); done; echo "\${out[*]}"; }

    {
        # Throwaway build file: a failed task is rerun, so no journal or fsync is needed
        echo "PRAGMA journal_mode = OFF;"
        echo "PRAGMA synchronous = OFF;"
        echo "PRAGMA temp_store = MEMORY;"
        echo "PRAGMA cache_size = -2000000;"
        echo ".bail on"
        echo "BEGIN;"
        # mgnifam by header name, so its schema column order is free and unlisted columns keep their DEFAULT
        c=\$(head -n 1 mgnifam.csv | tr -d '\\r')
        echo ".import --csv mgnifam.csv temp_mgnifam"
        echo "INSERT INTO mgnifam (\$c) SELECT \$(nullif "\$c") FROM temp_mgnifam;"
        echo "DROP TABLE temp_mgnifam;"
        # Child CSVs name their mgnifam_id column 'id', so these are imported by position
        for t in mgnifam_pfams mgnifam_funfams mgnifam_folds mgnifam_model_pfams; do
            [ -f "\$t.csv" ] || continue
            c=\$(cols "\$t")
            echo "CREATE TEMP TABLE temp_\$t (\$c);"
            echo ".import --csv --skip 1 \$t.csv temp_\$t"
            echo "INSERT INTO \$t (\$c) SELECT \$(nullif "\$c") FROM temp_\$t;"
            echo "DROP TABLE temp_\$t;"
        done
        echo "COMMIT;"
    } > import.sql

    sqlite3 "${prefix}.sqlite3" < import.sql > /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -L "${db}" "${prefix}.sqlite3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}
