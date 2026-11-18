-- One-off migration of an existing MGnifams sqlite DB: adds the families.tsv.gz provenance columns to mgnifam
-- (appended last; bin/load_families_tsv.py fills them by name). Guarded so reruns are no-ops:
--
--   db=mgnifams.sqlite3
--   [ "$(sqlite3 "$db" "SELECT count(*) FROM pragma_table_info('mgnifam') WHERE name = 'status'")" = 0 ] \
--       && sqlite3 -bail "$db" < assets/sql/migrate_schema_family_status.sql
--
-- Appended columns sit after the blobs, so reading them walks each row's overflow pages; a DB rebuilt from
-- assets/sql/db_schema.sqlite has them before the blobs.

BEGIN IMMEDIATE;

ALTER TABLE mgnifam ADD COLUMN status TEXT;
ALTER TABLE mgnifam ADD COLUMN not_updated_reason TEXT;
ALTER TABLE mgnifam ADD COLUMN first_release TEXT;
ALTER TABLE mgnifam ADD COLUMN model_release TEXT;
ALTER TABLE mgnifam ADD COLUMN model_proteins_release TEXT;
ALTER TABLE mgnifam ADD COLUMN members_release TEXT;
ALTER TABLE mgnifam ADD COLUMN members_proteins_release TEXT;

COMMIT;
