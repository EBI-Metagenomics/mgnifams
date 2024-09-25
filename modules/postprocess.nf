process QUERY_MGNPROTEIN_DB {
    publishDir "${params.outDir}", mode: "copy"
    container "quay.io/microbiome-informatics/mgnifams:2.0.0"

    input:
    path config_file
    path family_proteins_file

    output:
    path "post-processing"

    """
    python3 ${params.scriptDir}/post-processing/query_mgnprotein_db.py ${config_file} ${family_proteins_file}
    """
}

process PARSE_BIOMES {
    publishDir "${params.outDir}/post-processing", mode: "copy"

    input:
    path query_results

    output:
    path "biome_results"

    """
    python3 ${params.scriptDir}/post-processing/parse_biomes.py ${query_results}
    """
}

process CHUNK_QUERY_RESULTS_FOLDER {
    input:
    path post_processing
    val chunk_size

    output:
    path "chunk_folder*"

    """
    input_folder="${post_processing}/query_results"  # Path to the folder with the files
    chunk_number=${chunk_size}                        # Number of chunks you want
    output_base="chunk_folder"                        # Base name for the output folders (e.g., chunk_folder_1, chunk_folder_2, etc.)

    # Create output folders
    for ((i=1; i<=chunk_number; i++)); do
        mkdir -p "\${output_base}_\${i}"
    done

    # Find all files and symlinks in the input folder and split them into chunks
    file_count=\$(find "\$input_folder" \\( -type f -o -type l \\) | wc -l)
    files_per_chunk=\$(( (file_count + chunk_number - 1) / chunk_number ))  # Calculate files per chunk (ceil rounding)

    # Distribute files (and symlinks) into chunk folders
    i=1
    find "\$input_folder" \\( -type f -o -type l \\) | while read -r file; do
        # Determine the chunk folder to move the file/symlink to
        chunk=\$(( (i - 1) / files_per_chunk + 1 ))
        cp "\$file" "\${output_base}_\${chunk}/"
        ((i++))
    done
    """
}

process PARSE_DOMAINS {
    // publishDir "${params.outDir}/post-processing", mode: "copy"

    input:
    path chunk
    path query_results
    path refined_families

    output:
    path "domain_results"

    """
    python3 ${params.scriptDir}/post-processing/parse_domains.py ${chunk} ${query_results} ${refined_families}
    """
}

process INITIATE_SQLITE {
    container "docker://nouchka/sqlite3"

    input:
    path schema_file
    path tables
    
    output:
    path "mgnifams.sqlite3"

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

process APPEND_BLOBS {
    publishDir "${params.outDir}/post-processing", mode: "copy"
    label "venv"

    input:
    path db
    path cif_ch   , stageAs: "output/structures/cif/*"
    tuple val(meta2), path(seed_msa_ch, stageAs: "output/families/*")
    tuple val(meta3), path(msa_ch     , stageAs: "output/families/*")
    path hmm_ch   , stageAs: "output/families/hmm/*"
    path rf_ch    , stageAs: "output/families/rf/*"
    path biome_ch , stageAs: "output/post-processing/*"
    path domain_ch, stageAs: "output/post-processing/*"
    val update_blobs_from_row_id
    
    output:
    path "${db}"

    """
    python3 ${params.scriptDir}/post-processing/append_blobs_sqlite.py ${db} output families ${update_blobs_from_row_id}
    """
}
