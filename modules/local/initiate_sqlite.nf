process INITIATE_SQLITE {
    container "docker://nouchka/sqlite3"

    input:
    path schema_file
    path tables
    
    output:
    path "mgnifams.sqlite3"

    script:
    """
    sqlite3 mgnifams.sqlite3 < ${schema_file}

    # Import mgnifam.csv directly with NULLs for missing columns
    sqlite3 mgnifams.sqlite3 <<EOF
    .mode csv
    .import '${tables}/mgnifam.csv' mgnifam
    .exit
    EOF

    # Create temporary tables for the other CSV files
    sqlite3 mgnifams.sqlite3 <<EOF
    CREATE TEMP TABLE temp_mgnifam_proteins (
        mgnifam_id INTEGER,
        protein INTEGER,
        region TEXT
    );

    CREATE TEMP TABLE temp_mgnifam_pfams (
        mgnifam_id INTEGER,
        rank INTEGER,
        pfam_id VARCHAR(8),
        pfam_hit TEXT,
        query_hmm_range TEXT,
        template_hmm_range TEXT,
        e_value REAL
    );

    CREATE TEMP TABLE temp_mgnifam_folds (
        mgnifam_id INTEGER,
        target_structure TEXT,
        aligned_length INTEGER,
        query_start INTEGER,
        query_end INTEGER,
        target_start INTEGER,
        target_end INTEGER,
        e_value REAL
    );
    .exit
    EOF

    # Import data into temporary tables
    sqlite3 mgnifams.sqlite3 <<EOF
    .mode csv
    .import '${tables}/mgnifam_proteins.csv' temp_mgnifam_proteins
    .import '${tables}/mgnifam_pfams.csv' temp_mgnifam_pfams
    .import '${tables}/mgnifam_folds.csv' temp_mgnifam_folds
    .exit
    EOF

    # Insert data from temporary tables into the main tables
    sqlite3 mgnifams.sqlite3 <<EOF
    INSERT INTO mgnifam_proteins (mgnifam_id, protein, region)
    SELECT mgnifam_id, protein, region FROM temp_mgnifam_proteins;

    INSERT INTO mgnifam_pfams (mgnifam_id, rank, pfam_id, pfam_hit, query_hmm_range, template_hmm_range, e_value)
    SELECT mgnifam_id, rank, pfam_id, pfam_hit, query_hmm_range, template_hmm_range, e_value FROM temp_mgnifam_pfams;

    INSERT INTO mgnifam_folds (mgnifam_id, target_structure, aligned_length, query_start, query_end, target_start, target_end, e_value)
    SELECT mgnifam_id, target_structure, aligned_length, query_start, query_end, target_start, target_end, e_value FROM temp_mgnifam_folds;

    DROP TABLE temp_mgnifam_proteins;
    DROP TABLE temp_mgnifam_pfams;
    DROP TABLE temp_mgnifam_folds;
    .exit
    EOF
    """
}
