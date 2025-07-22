process IMPORT_QUERIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d41320ac1ca5f0a626982296d23ded50376a966b9e8240aa50dba2014a805bf5/data':
        'community.wave.seqera.io/library/sqlite:3.48.0--48957425ca78aa09' }"

    input:
    tuple val(meta) , path(pipeline_results)
    tuple val(meta2), path(db)
    
    output:
    tuple val(meta), path(db), emit: db
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Import mgnifam.csv directly with NULLs for missing columns
    sqlite3 ${db} <<EOF
    .mode csv
    .import 'mgnifam.csv' mgnifam
    .exit
    EOF

    # Create temporary tables for the other CSV files
    sqlite3 ${db} <<EOF
    CREATE TEMP TABLE temp_mgnifam_funfams (
        mgnifam_id INTEGER,
        funfam TEXT,
        e_value REAL,
        score REAL,
        hmm_from INTEGER,
        hmm_to INTEGER,
        ali_from INTEGER,
        ali_to INTEGER,
        env_from INTEGER,
        env_to INTEGER,
        acc REAL
    );

    CREATE TEMP TABLE temp_mgnifam_folds (
        mgnifam_id INTEGER,
        fold TEXT,
        aligned_length INTEGER,
        q_start INTEGER,
        q_end INTEGER,
        t_start INTEGER,
        t_end INTEGER,
        e_value REAL
    );

    CREATE TEMP TABLE temp_mgnifam_pfams (
        mgnifam_id INTEGER,
        pfam TEXT,
        name TEXT,
        description TEXT,
        prob REAL,
        e_value REAL,
        length INTEGER,
        query_hmm TEXT,
        template_hmm TEXT
    );
    .exit
    EOF

    # Import data into temporary tables
    sqlite3 ${db} <<EOF
    .mode csv
    .import 'mgnifam_funfams.csv' temp_mgnifam_funfams
    .import 'mgnifam_folds.csv' temp_mgnifam_folds
    .import 'mgnifam_pfams.csv' temp_mgnifam_pfams
    .exit
    EOF

    # Insert data from temporary tables into the main tables
    sqlite3 ${db} <<EOF
    INSERT INTO mgnifam_funfams (mgnifam_id, funfam, e_value, score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc)
    SELECT id, funfam, e_value, score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc FROM temp_mgnifam_funfams;

    INSERT INTO mgnifam_folds (mgnifam_id, fold, aligned_length, q_start, q_end, t_start, t_end, e_value)
    SELECT id, fold, aligned_length, q_start, q_end, t_start, t_end, e_value FROM temp_mgnifam_folds;

    INSERT INTO mgnifam_pfams (mgnifam_id, pfam, name, description, prob, e_value, length, query_hmm, template_hmm)
    SELECT id, pfam, name, description, prob, e_value, length, query_hmm, template_hmm FROM temp_mgnifam_pfams;

    DROP TABLE temp_mgnifam_funfams;
    DROP TABLE temp_mgnifam_folds;
    DROP TABLE temp_mgnifam_pfams;
    .exit
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}
