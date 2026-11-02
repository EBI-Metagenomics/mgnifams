# AGENTS.md

Guidance for coding agents (Claude Code, Codex and others) working in this repository.

## Overview

MGnifams is a Nextflow (DSL2) pipeline that converts metagenomics-derived amino acid sequences into protein families. It follows nf-core standards and is deployed at EMBL-EBI for the MGnify platform.

`dev` is at `3.0.0`, the release candidate: the version was bumped before the SLURM update run because `workflow.manifest.version` is written into `update_info.csv` and from there into `release_history.pipeline_version` of MGnifams 1.1. After the tag, `dev` goes to `3.1.0dev`. The release is a major because it holds breaking changes (see the CHANGELOG). `CHANGELOG.md` entries link their PR, newest first, with same-PR points of one category grouped as a sub-list; tool version changes go in the `Dependencies` table.

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

2. **GENERATE_NONREDUNDANT_FAMILIES** (`subworkflows/local/generate_nonredundant_families/`) — Core algorithm in `bin/generate_families.py` (to be replaced by nf-core `mgnifam/generatefamilies` 4.0.0, [#68](https://github.com/EBI-Metagenomics/mgnifams/issues/68)). Iteratively builds families: seed MSA (pyfamsa) → HMM (pyhmmer) → recruit sequences → align (pyhmmer/hmmalign) → trim (pytrimal), repeated up to 3 times. Redundancy removed via `subworkflows/local/remove_redundancy/` using `bin/identify_redundant_fams.py`.

3. **PREDICT_STRUCTURES** (`subworkflows/local/predict_structures/`) — ESMFold protein structure prediction with GPU support. CUDA OOM failures are caught and re-run on CPU (`bin/extract_cuda_failed.py`). The CPU re-run (`RUN_ESMFOLD_CPU`) is an alias of `RUN_ESMFOLD`, so it carries the static `process_gpu` label; `conf/base.config` clears its `accelerator` (a config `withName` cannot remove a label). Site configs must therefore request GPUs from `task.accelerator` (`clusterOptions = { task.accelerator ? '--gres=gpu:1' : '' }`), because the Nextflow SLURM executor ignores `accelerator` itself. Outputs PDB/CIF files with pLDDT and pTM scores.

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
exactly one row: `mgnify_proteins_sequences`, `mgnify_proteins_clusters` and `mgnify_proteins_pfam` (MGnify proteins parquet files),
`mgnifams_hmms` (HMM library with numeric `NAME`s), an optional `mgnprotein_db_config` and `mgnifams_families` (the
previous release's FTP `families.tsv.gz`).

1. **UPDATE_FAMILIES** (`subworkflows/local/update_families/`) — `EXTRACT_UNANNOTATED_PARQUET_SLICES` slices known Pfam
   domains off the MGnify90 cluster representatives only (`cluster_rep` of `mgy_clusters.parquet`, as in 1.0; parquet row
   groups split into `--parquet_chunks`), `SPLIT_HMM_LIB` chunks the library by
   `--hmm_chunk_size`, nf-core `mgnifam/updatefamilies` (`--skip_refine`) recruits and aligns, and
   `POOL_UPDATED_FAMILIES` pools the chunks (`update_families/updated_delta.csv` holds the outcome per family), failing
   if a chunk's CSV header is not the expected one. `family_metadata.csv` keeps the mgnifam 4.0.0 header
   (`family_id,converged,seed_msa_size,full_msa_size,rep_protein,rep_region,rep_length,consensus_length,rep_sequence,consensus_sequence`;
   `converged`/`seed_msa_size` empty under `--skip_refine`). Default mode still writes the old 8-column header; readers
   (`bin/export_mgnifams.py`, `bin/export_families_tsv.py`) select columns by name.
   `EXPORT_FAMILIES_TSV` then writes the release's `update_families/families.tsv.gz` from `mgnifams_families`, the delta
   and `family_metadata.csv` (`bin/export_families_tsv.py`; fails unless the previous and updated family sets match).
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

- **`nextflow.config`** — Main config with 80+ parameters, profiles (docker/singularity/conda/slurm/gpu/test), and the nf-schema 2.7.2 plugin
- **`conf/base.config`** — Resource tiers: single (1 CPU/6GB/4h), low, medium, high, plus custom GENERATE_FAMILIES (4 CPU, 1400 GB/72 h × attempt, 3 retries)
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

## FTP releases

`ftp/` is the draft of the MGnifams FTP layout (`ftp/README.md` documents it): only READMEs and small metadata files are
committed. `update_mgnifams` writes the release's `families.tsv.gz` (`EXPORT_FAMILIES_TSV`, columns in `ftp/README.md`).
**TODO (with [#68](https://github.com/EBI-Metagenomics/mgnifams/issues/68)):** `run_mgnifams_pipeline` (`workflows/mgnifams.nf`)
does not yet; 1.0's was a one-off export from the prod DB. Generate it there (all `new`, `first_release` = `model_release` =
`members_release`) when that workflow switches to `mgnifam/generatefamilies`.

## Linting & hooks

- `prek install` once, then `prek run --all-files` (config: `.pre-commit-config.yaml`): ruff check/format for `bin/`,
  prettier for YAML/JSON/Markdown, and `nextflow lint`.
- `nextflow lint . -o concise | grep -E '^(Warn|Error)' | grep -v nf-core/` must print nothing: 0 warnings in
  pipeline-owned files (`main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/`).
- A workflow with one emit must leave it unnamed. A bare identifier (`emit: ch_x`) still counts as named, so emit an
  expression instead (e.g. `hh_mode == "hhblits" ? HHSUITE_HHBLITS.out.hhr : HHSUITE_HHSEARCH.out.hhr`) and read it
  with `.out` in the caller.

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

- **Software versions:** nf-core modules emit `[process, tool, version]` tuples to the `versions` topic; local modules
  still emit `versions.yml` into `ch_versions`. Never mix topic outputs into `ch_versions`: each top-level workflow
  (`workflows/mgnifams.nf`, `workflows/update_mgnifams.nf`) reads `channel.topic("versions")` once, before MultiQC,
  and merges it with `softwareVersionsToYAML(ch_versions)`.
- **MultiQC** (nf-core module) takes one tuple `[meta, files, configs, logo, replace_names, sample_names]`; the custom
  `--multiqc_config` goes after `assets/multiqc_config.yml` so it takes precedence.
- **Patches:** only `s4pred/runmodel` is patched (`s4pred-runmodel.diff`, FASTA header cleanup). After
  `nf-core modules update`, re-apply and regenerate it with `nf-core modules patch`. Lines starting with `---` inside a
  patch break nf-core tools' parser.
- `conf/containers_*.config` are written by nf-core tools (`nf-core modules patch/update`); they are committed but not
  included by `nextflow.config`.
- Local module/subworkflow tests pull test data from nf-core/test-datasets, branch `proteinfamilies`
  (`params.pipelines_testdata_base_path`), e.g. `test_data/mgnifams_input_small.faa`.
