-- Merges an update_mgnifams delta DB (db/<sample>_update.sqlite3) into a production MGnifams DB, in one transaction.
-- Prod must be migrated first (assets/migrate_schema_seed_size_tmscores.sql). Run with the delta attached:
--
--   sqlite3 -bail prod.sqlite3 -cmd "ATTACH 'delta.sqlite3' AS delta" < assets/merge_update_delta.sql
--
-- Any failed check or statement exits before COMMIT (-bail), so prod is left untouched.
-- Kept from prod: consensus, hmm_length, converged, seed_msa_blob, hmm_blob, rf_blob, seed_size, mgnifam_model_pfams,
-- and the TM columns / biome_blob unless the delta's update_info says they were computed.

PRAGMA foreign_keys = ON;

BEGIN IMMEDIATE;

-- A failed check aborts with "CHECK constraint failed: <constraint name>"
CREATE TEMP TABLE merge_check (
    delta_complete INTEGER DEFAULT 1 CONSTRAINT delta_lacks_update_info_child_tables CHECK (delta_complete IS 1),
    ids_in_prod    INTEGER DEFAULT 1 CONSTRAINT delta_has_ids_missing_from_prod CHECK (ids_in_prod IS 1),
    foreign_keys   INTEGER DEFAULT 1 CONSTRAINT merge_breaks_foreign_keys CHECK (foreign_keys IS 1)
);

-- Preflight: a complete delta (update_info written last by UPDATE_SQLITE_BLOBS_STAGED) of existing families only
INSERT INTO temp.merge_check (delta_complete) SELECT
    (SELECT value FROM delta.update_info WHERE key = 'child_tables') = 'pfams,funfams,folds';
INSERT INTO temp.merge_check (ids_in_prod) SELECT
    NOT EXISTS (SELECT 1 FROM delta.mgnifam WHERE id NOT IN (SELECT id FROM main.mgnifam));

UPDATE main.mgnifam AS m SET
    full_size      = d.full_size,
    protein_rep    = d.protein_rep,
    rep_region     = d.rep_region,
    rep_length     = d.rep_length,
    rep_sequence   = d.rep_sequence,
    plddt          = d.plddt,
    ptm            = d.ptm,
    helix_percent  = d.helix_percent,
    strand_percent = d.strand_percent,
    coil_percent   = d.coil_percent,
    cif_blob       = d.cif_blob,
    domain_blob    = d.domain_blob,
    s4pred_blob    = d.s4pred_blob
FROM delta.mgnifam AS d
WHERE m.id = d.id;

UPDATE main.mgnifam AS m SET
    inside_percent         = d.inside_percent,
    membrane_alpha_percent = d.membrane_alpha_percent,
    outside_percent        = d.outside_percent,
    signal_percent         = d.signal_percent,
    membrane_beta_percent  = d.membrane_beta_percent,
    periplasm_percent      = d.periplasm_percent,
    tm_blob                = d.tm_blob
FROM delta.mgnifam AS d
WHERE m.id = d.id AND (SELECT value FROM delta.update_info WHERE key = 'tm_computed') = 'true';

UPDATE main.mgnifam AS m SET
    biome_blob = d.biome_blob
FROM delta.mgnifam AS d
WHERE m.id = d.id AND (SELECT value FROM delta.update_info WHERE key = 'biome_computed') = 'true';

-- Child tables are replaced for the updated families; the AUTOINCREMENT id is left to prod
DELETE FROM main.mgnifam_pfams   WHERE mgnifam_id IN (SELECT id FROM delta.mgnifam);
DELETE FROM main.mgnifam_funfams WHERE mgnifam_id IN (SELECT id FROM delta.mgnifam);
DELETE FROM main.mgnifam_folds   WHERE mgnifam_id IN (SELECT id FROM delta.mgnifam);

INSERT INTO main.mgnifam_pfams (mgnifam_id, pfam, name, e_value, score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc)
SELECT mgnifam_id, pfam, name, e_value, score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc FROM delta.mgnifam_pfams;

INSERT INTO main.mgnifam_funfams (mgnifam_id, funfam, e_value, score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc)
SELECT mgnifam_id, funfam, e_value, score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, acc FROM delta.mgnifam_funfams;

INSERT INTO main.mgnifam_folds (mgnifam_id, fold, aligned_length, q_start, q_end, t_start, t_end, e_value, aln_tmscore, q_tmscore, t_tmscore)
SELECT mgnifam_id, fold, aligned_length, q_start, q_end, t_start, t_end, e_value, aln_tmscore, q_tmscore, t_tmscore FROM delta.mgnifam_folds;

INSERT INTO temp.merge_check (foreign_keys) SELECT NOT EXISTS (SELECT 1 FROM main.pragma_foreign_key_check);

COMMIT;
