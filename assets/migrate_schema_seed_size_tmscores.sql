-- One-off migration of an existing MGnifams sqlite DB: adds mgnifam.seed_size and the Foldseek TM-score
-- columns of mgnifam_folds (appended last, like the fresh schema in assets/data/db_schema.sqlite).
-- seed_size is backfilled from seed_msa_blob (aligned FASTA); old fold rows keep NULL TM-scores until
-- their families are updated. Run once, before the first update merge; it is guarded so reruns are no-ops:
--
--   db=mgnifams.sqlite3
--   [ "$(sqlite3 "$db" "SELECT count(*) FROM pragma_table_info('mgnifam') WHERE name = 'seed_size'")" = 0 ] \
--       && sqlite3 -bail "$db" < assets/migrate_schema_seed_size_tmscores.sql
--
-- Without the guard, a second run fails on "duplicate column name" and -bail rolls it back.

BEGIN IMMEDIATE;

ALTER TABLE mgnifam ADD COLUMN seed_size INTEGER;
ALTER TABLE mgnifam_folds ADD COLUMN aln_tmscore REAL;
ALTER TABLE mgnifam_folds ADD COLUMN q_tmscore REAL;
ALTER TABLE mgnifam_folds ADD COLUMN t_tmscore REAL;

UPDATE mgnifam
SET seed_size = length(CAST(seed_msa_blob AS TEXT)) - length(replace(CAST(seed_msa_blob AS TEXT), '>', ''))
WHERE seed_msa_blob IS NOT NULL;

COMMIT;
