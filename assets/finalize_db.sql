-- Search flags, indexes and planner stats of a full MGnifams DB. Idempotent: FINALIZE_SQLITE runs it after the import
-- (init_mgnifams_db), and it can be rerun by hand on any DB with the has_* columns.
-- has_structure = has a Foldseek hit (mgnifam_folds); every family has a structure in cif_blob.

BEGIN;

-- Child indexes first: the flag subqueries below use them
CREATE INDEX IF NOT EXISTS idx_pfams_mgnifam          ON mgnifam_pfams (mgnifam_id);
CREATE INDEX IF NOT EXISTS idx_funfams_mgnifam        ON mgnifam_funfams (mgnifam_id);
CREATE INDEX IF NOT EXISTS idx_folds_mgnifam          ON mgnifam_folds (mgnifam_id);
CREATE INDEX IF NOT EXISTS idx_model_pfams_mgnifam    ON mgnifam_model_pfams (mgnifam_id);

UPDATE mgnifam SET
    has_pfam       = EXISTS (SELECT 1 FROM mgnifam_pfams       WHERE mgnifam_id = mgnifam.id),
    has_funfam     = EXISTS (SELECT 1 FROM mgnifam_funfams     WHERE mgnifam_id = mgnifam.id),
    has_model_pfam = EXISTS (SELECT 1 FROM mgnifam_model_pfams WHERE mgnifam_id = mgnifam.id),
    has_structure  = EXISTS (SELECT 1 FROM mgnifam_folds       WHERE mgnifam_id = mgnifam.id);

CREATE INDEX IF NOT EXISTS idx_mgnifam_full_size      ON mgnifam (full_size);
CREATE INDEX IF NOT EXISTS idx_mgnifam_seed_size      ON mgnifam (seed_size);
CREATE INDEX IF NOT EXISTS idx_mgnifam_rep_length     ON mgnifam (rep_length);
CREATE INDEX IF NOT EXISTS idx_mgnifam_hmm_length     ON mgnifam (hmm_length);
CREATE INDEX IF NOT EXISTS idx_mgnifam_plddt          ON mgnifam (plddt);
CREATE INDEX IF NOT EXISTS idx_mgnifam_ptm            ON mgnifam (ptm);
CREATE INDEX IF NOT EXISTS idx_mgnifam_helix          ON mgnifam (helix_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_strand         ON mgnifam (strand_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_coil           ON mgnifam (coil_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_inside         ON mgnifam (inside_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_membrane_alpha ON mgnifam (membrane_alpha_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_outside        ON mgnifam (outside_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_signal         ON mgnifam (signal_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_membrane_beta  ON mgnifam (membrane_beta_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_periplasm      ON mgnifam (periplasm_percent);
CREATE INDEX IF NOT EXISTS idx_mgnifam_has_pfam       ON mgnifam (has_pfam);
CREATE INDEX IF NOT EXISTS idx_mgnifam_has_funfam     ON mgnifam (has_funfam);
CREATE INDEX IF NOT EXISTS idx_mgnifam_has_model_pfam ON mgnifam (has_model_pfam);
CREATE INDEX IF NOT EXISTS idx_mgnifam_has_structure  ON mgnifam (has_structure);

COMMIT;

ANALYZE;
