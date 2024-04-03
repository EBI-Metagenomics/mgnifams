CREATE TABLE mgnifam (
    id INT PRIMARY KEY,
    family_size INT,
    protein_rep INT,
    rep_region TEXT,
    converged BOOLEAN,
    cif_file TEXT,
    seed_msa_file TEXT,
    msa_file TEXT,
    hmm_file TEXT,
    biomes_file TEXT,
    domain_architecture_file TEXT
);

CREATE TABLE mgnifam_proteins (
    id SERIAL PRIMARY KEY,
    mgnifam_id INT REFERENCES mgnifam(id),
    protein INT,
    region TEXT
);

CREATE TABLE mgnifam_pfams (
    id SERIAL PRIMARY KEY,
    mgnifam_id INT REFERENCES mgnifam(id),
    rank INT,
    pfam_id VARCHAR(8),
    pfam_hit TEXT,
    query_hmm_range TEXT,
    template_hmm_range TEXT,
    e_value DOUBLE PRECISION
);

CREATE TABLE mgnifam_folds (
    id SERIAL PRIMARY KEY,
    mgnifam_id INT REFERENCES mgnifam(id),
    target_structure TEXT,
    aligned_length INT,
    query_start INT,
    query_end INT,
    target_start INT,
    target_end INT,
    e_value DOUBLE PRECISION
);
