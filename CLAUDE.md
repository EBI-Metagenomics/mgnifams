# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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
```

## Testing

```bash
# End-to-end nf-test
nf-test test tests/default.nf.test --profile +singularity,test
```

The test profile requires paths to external databases set in your local config:
- `esmfold_db`, `esmfold_params_path` — ESMFold model weights
- `pfam_path` — Pfam-A HMM database (gzipped)
- `funfams_path` — FunFams HMM library (gzipped)
- `hhdb_path` — HHsuite Pfam database directory
- `foldseek_db_path` — Foldseek database directory

## Pipeline Architecture

Three modes controlled by `--mode` parameter:

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
Initializes SQLite schema (`assets/data/db_schema.sqlite`) and imports pipeline results.

### `update_mgnifams_db` — `workflows/update_db.nf`
Queries MGnify proteins PostgreSQL DB (`bin/query_mgnprotein_db.py`) to enrich families with biome and domain architecture data, then updates the SQLite via `bin/update_sqlite_blobs.py`.

## Key Configuration

- **`nextflow.config`** — Main config with 80+ parameters, profiles (docker/singularity/conda/slurm/gpu/test), and nf-schema v2.3.0 plugin
- **`conf/base.config`** — Resource tiers: single (1 CPU/6GB/4h), low, medium, high, plus custom GENERATE_FAMILIES (4 CPU/400GB/35h, 5 retries)
- **`conf/modules.config`** — Per-process resource overrides and tool parameters
- **`nextflow_schema.json`** — Full parameter validation schema
- **`assets/schema_input.json`** — Samplesheet validation

The `conf/slurm.config` and `conf/local.config` files are gitignored — create them locally.

## Code Layout

- **`bin/`** — Python helper scripts called by Nextflow processes. `generate_families.py` is the core algorithm (~26KB).
- **`modules/local/`** — Custom Nextflow process definitions
- **`modules/nf-core/`** — Standard nf-core modules (don't modify directly)
- **`subworkflows/local/`** — Custom subworkflow logic
- **`subworkflows/nf-core/`** — Standard nf-core subworkflows
- **`workflows/`** — Top-level workflow entry points

## nf-core Conventions

This pipeline follows nf-core DSL2 conventions. When adding modules, use `nf-core modules install` or follow patterns in `modules/local/`. The `modules.json` tracks nf-core module versions. nf-test is used for testing; nf-core module tests in `modules/nf-core/**/tests/` are excluded from the local nf-test config.
