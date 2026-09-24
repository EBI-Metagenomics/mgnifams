# AGENTS.md

Guidance for coding agents (Claude Code, Codex and others) working in this repository.

## Overview

MGnifams is a Nextflow (DSL2) pipeline that converts metagenomics-derived amino acid sequences into protein families. It follows nf-core standards and is deployed at EMBL-EBI for the MGnify platform.

## Running the Pipeline

```bash
# Local execution (test profile)
nextflow run . -c ../conf/local.config --input input/samplesheet_test.csv -profile test,local,singularity -resume

# HPC/SLURM with GPU
nextflow run . -c ../conf/slurm.config --input input/samplesheet_test.csv -profile test,slurm,singularity,gpu -resume -with-tower

# Initialize SQLite database
nextflow run . -c ../conf/slurm.config --input input/samplesheet_init_db.csv --mode init_mgnifams_db --outdir '/path/to/output_db' -profile slurm,singularity -resume

# Update database from MGnify proteins DB
nextflow run . -c ../conf/slurm.config --input input/samplesheet_update_db.csv --mode update_mgnifams_db --outdir '/path/to/output_db' -profile slurm,singularity -resume

# Update existing families with new MGnify proteins
nextflow run . -c ../conf/local.config -profile test_update_mgnifams,local,singularity -resume

# Delta DB from an update_mgnifams outdir (default: the tests/data/update_results fixture; run from the repo root)
nextflow run . -profile test_post_update_mgnifams_update_db,singularity
```

## Testing

```bash
# End-to-end nf-test
nf-test test tests/default.nf.test --profile +singularity,test
nf-test test tests/update_mgnifams.nf.test --profile +singularity
nf-test test tests/post_update_mgnifams_update_db.nf.test --profile +singularity

# Module tests and script self-checks
nf-test test modules/local/<module> --profile +singularity
python3 bin/test_<script>.py
```

`tests/update_mgnifams.nf.test` has a self-contained `-stub` variant. nf-test cannot load `../conf/local.config`, so for its real test (as for `tests/default.nf.test`) add the database paths to `conf/test_update_mgnifams.config` locally, without committing them.
The real pipeline tests run ESMFold, which needs most of a 30 GB machine on CPU.

The test profile requires paths to external databases set in your local config:

- `esmfold_db`, `esmfold_params_path` — ESMFold model weights
- `pfam_path` — Pfam-A HMM database (gzipped)
- `funfams_path` — FunFams HMM library (gzipped)
- `hhdb_path` — HHsuite Pfam database directory
- `foldseek_db_path` — Foldseek database directory

## Pipeline Architecture

Five modes controlled by `--mode` parameter:

### `run_mgnifams_pipeline` (default) — `workflows/mgnifams.nf`

Five sequential subworkflows:

1. **SETUP_CLUSTERS** (`subworkflows/local/setup_clusters/`) — Extract unannotated sequences from MGnify CSV (or use FASTA directly via `--fasta_input_mode`), filter by length, quality check with seqkit, cluster with mmseqs linclust, distribute into chunks.

2. **GENERATE_NONREDUNDANT_FAMILIES** (`subworkflows/local/generate_nonredundant_families/`) — Core algorithm in `bin/generate_families.py`. Iteratively builds families: seed MSA (pyfamsa) → HMM (pyhmmer) → recruit sequences → align (pyhmmer/hmmalign) → trim (pytrimal), repeated up to 3 times. Redundancy removed via `subworkflows/local/remove_redundancy/` using `bin/identify_redundant_fams.py`.

3. **PREDICT_STRUCTURES** (`subworkflows/local/predict_structures/`) — ESMFold protein structure prediction with GPU support. CUDA OOM failures are caught and re-run on CPU (`bin/extract_cuda_failed.py`). Outputs PDB/CIF files with pLDDT and pTM scores.

4. **ANNOTATE_FAMILIES** (`subworkflows/local/annotate_families/`) — Three parallel annotation streams:
   - `ANNOTATE_REPS`: secondary structure (s4pred), transmembrane (deeptmhmm), Pfam/FunFams hmmsearch
   - `ANNOTATE_MODELS`: family-level Pfam via hhsuite/hhsearch
   - `ANNOTATE_STRUCTURES`: structural homologs via foldseek against PDB/AlphaFoldDB

5. **EXPORT_DATA** (`subworkflows/local/export_data/`) — Tabular output for the web interface and MultiQC reports.

### `init_mgnifams_db` — `workflows/init_db.nf`

Initializes SQLite schema (`assets/data/db_schema.sqlite`), imports pipeline results, then `FINALIZE_SQLITE` runs
`assets/finalize_db.sql` (`has_*` flags, indexes, `ANALYZE`).

### `update_mgnifams_db` — `workflows/update_db.nf`

Queries MGnify proteins PostgreSQL DB (`bin/query_mgnprotein_db.py`) to enrich families with biome and domain architecture data, then updates the SQLite via `bin/update_sqlite_blobs.py`.

### `update_mgnifams` — `workflows/update_mgnifams.nf`

Refreshes the full MSAs and representatives of existing families; seed MSAs and HMMs do not change. The samplesheet has
exactly one row: `mgnify_proteins_sequences` and `mgnify_proteins_pfam` (MGnify proteins parquet files),
`mgnifams_hmms` (HMM library with numeric `NAME`s) and an optional `mgnprotein_db_config`.

1. **UPDATE_FAMILIES** (`subworkflows/local/update_families/`) — `EXTRACT_UNANNOTATED_PARQUET_SLICES` slices known Pfam
   domains off the proteins (parquet row groups split into `--parquet_chunks`), `SPLIT_HMM_LIB` chunks the library by
   `--hmm_chunk_size`, nf-core `mgnifam/updatefamilies` (`--skip_refine`) recruits and aligns, and
   `POOL_UPDATED_FAMILIES` pools the chunks (`update_families/updated_delta.csv` holds the outcome per family).
2. On the successful families: `PREDICT_STRUCTURES`, `ANNOTATE_REPS`, `ANNOTATE_STRUCTURES`, `EXPORT_DATA`; domain
   architectures from the Pfam parquet (`BUILD_PARQUET_DOMAIN_QUERIES` → `PARSE_DOMAINS`); biomes only with
   `mgnprotein_db_config`. `--run_alphafold2` adds ColabFold predictions from the full MSAs (`structures/alphafold2/`, not used downstream).
3. `update_families/update_info.csv` records `tm_computed` / `biome_computed` for the DB mode below. No sqlite work.

### `post_update_mgnifams_update_db` — `workflows/post_update_mgnifams_update_db.nf`

Delta DB `db/<sample>_update.sqlite3` from an `update_mgnifams` outdir (samplesheet `sample,results_folder`, one row):
`INIT_SQLITE` → `IMPORT_QUERIES` → `UPDATE_SQLITE_BLOBS_STAGED` (fails instead of publishing a partial DB), staging the
published CSVs, CIFs, feature JSONs, domain/biome results and `family_ids.fasta` ids. Its `update_info` table comes from
`update_info.csv`. The pipeline never touches prod: `assets/merge_update_delta.sql` merges the delta (see the README
for which columns it overwrites and keeps). Test: `tests/post_update_mgnifams_update_db.nf.test` on the two-family
fixture `tests/data/update_results` (no external DBs, cheap).

## Key Configuration

- **`nextflow.config`** — Main config with 80+ parameters, profiles (docker/singularity/conda/slurm/gpu/test), and nf-schema v2.3.0 plugin
- **`conf/base.config`** — Resource tiers: single (1 CPU/6GB/4h), low, medium, high, plus custom GENERATE_FAMILIES (4 CPU/400GB/35h, 5 retries)
- **`conf/modules.config`** — Per-process resource overrides and tool parameters
- **`nextflow_schema.json`** — Full parameter validation schema
- **`assets/schema_input.json`** — Samplesheet validation

The `conf/slurm.config` and `conf/local.config` files are gitignored — create them locally.

## DB schema

`assets/data/db_schema.sqlite` (SQL text) now has `mgnifam.seed_size` and the Foldseek TM-scores
`mgnifam_folds.aln_tmscore/q_tmscore/t_tmscore`. `IMPORT_QUERIES` imports `mgnifam.csv` by header name (every
`export_mgnifams.py` header column must exist in the schema, checked by `bin/test_export_mgnifams.py`); columns missing
from the CSV keep their DEFAULT. Child CSVs are imported by position (their `id` column is the `mgnifam_id`).
Migrated prod DBs have `seed_size`/`hmm_length` last instead; that is fine, since `merge_update_delta.sql` uses column names.
The `has_*` flags are derived from the child tables: `finalize_db.sql` fills them (init DB), the merge recomputes them
for the delta families; the delta DB leaves them at 0. Fast PRAGMAs (journal/fsync off) only on throwaway build files, never prod.
`mgnifam.hmm_length` = `length(consensus)` (one consensus residue per match state), stored so the website can index/sort on it.
Existing DBs are migrated once with `assets/migrate_schema_seed_size_tmscores.sql`.
`release_history` has one row per MGnifams release in the DB (`release` is `MAJOR.MINOR`, independent of the MGnify Proteins
`YYYY_MM` release it records). `merge_update_delta.sql` creates it if missing and adds the update row from the delta's
`update_info` (`mgnifams_release`, `mgnify_proteins_release` come from the `update_mgnifams` params of the same names).
Discarded families are never in the delta, so the merge leaves them at their previous version.

## Linting & hooks

- `prek install` once, then `prek run --all-files` (config: `.pre-commit-config.yaml`): ruff check/format for `bin/`,
  prettier for YAML/JSON/Markdown, and `nextflow lint`.
- `nextflow lint . -o concise | grep -E '^(Warn|Error)' | grep -v nf-core/` must print nothing: 0 warnings in
  pipeline-owned files (`main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/`).

## Code Layout

- **`bin/`** — Python helper scripts called by Nextflow processes. `generate_families.py` is the core algorithm (~26KB).
- **`bin/test_*.py`** — Assert-based self-checks of the scripts next to them (`python3 bin/test_<name>.py`)
- **`modules/local/`** — Custom Nextflow process definitions
- **`assets/*.sql`** — Out-of-pipeline SQL for production DBs: schema migration, update-delta merge
- **`modules/nf-core/`** — Standard nf-core modules (don't modify directly)
- **`subworkflows/local/`** — Custom subworkflow logic
- **`subworkflows/nf-core/`** — Standard nf-core subworkflows
- **`workflows/`** — Top-level workflow entry points

## nf-core Conventions

This pipeline follows nf-core DSL2 conventions. When adding modules, use `nf-core modules install` or follow patterns in `modules/local/`. The `modules.json` tracks nf-core module versions. nf-test is used for testing; nf-core module tests in `modules/nf-core/**/tests/` are excluded from the local nf-test config.
