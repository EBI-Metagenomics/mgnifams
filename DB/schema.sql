CREATE TABLE mgnifam (
    id INT PRIMARY KEY,
    family_size INT,
    protein_rep INT,
    mask VARCHAR,
    cif_file TEXT,
    msa_file TEXT,
    seed_msa_file TEXT,
    hmm_file TEXT
);

CREATE TABLE mgnifam_pfams (
    id SERIAL PRIMARY KEY,
    mgnifam_id INT REFERENCES mgnifam(id),
    rank INT,
    pfam_hit VARCHAR,
    query_hmm_range VARCHAR,
    template_hmm_range VARCHAR,
    e_value DOUBLE PRECISION
);

CREATE TABLE mgnifam_folds (
    id SERIAL PRIMARY KEY,
    mgnifam_id INT REFERENCES mgnifam(id),
    target_structure VARCHAR,
    aligned_length INT,
    query_start INT,
    query_end INT,
    target_start INT,
    target_end INT,
    e_value DOUBLE PRECISION
);