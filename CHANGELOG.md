# ebi-metagenomics/mgnifams: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v3.0.0 - [unreleased]

### `Added`

- [#66](https://github.com/EBI-Metagenomics/mgnifams/pull/66) - `update_mgnifams` searches only the MGnify90 cluster representatives, as MGnifams 1.0 did, read from the new required samplesheet column `mgnify_proteins_clusters` (`mgy_clusters.parquet`).
- [#65](https://github.com/EBI-Metagenomics/mgnifams/pull/65) - Release provenance:
  - New `release_history` table: one row per MGnifams release held by the DB (release, date, type, MGnify Proteins release, pipeline version, family counts). `assets/merge_update_delta.sql` creates it if missing and records the update release; merging the same release twice fails.
  - New `update_mgnifams` parameters `--mgnifams_release` (default `1.1`) and `--mgnify_proteins_release` (default `2026_07`), validated as `MAJOR.MINOR` and `YYYY_MM` and written to `update_info.csv`.
  - `ftp/`: draft of the versioned MGnifams FTP layout (`releases/<MAJOR.MINOR>/`, `RELEASES.tsv`, per-release READMEs and `release_manifest.json`), with 1.0 filled in from the existing release.
- [#64](https://github.com/EBI-Metagenomics/mgnifams/pull/64) - Update mechanism and DB schema additions:
  - New `--mode update_mgnifams`: refreshes existing families against new MGnify proteins, read from the MGnify proteins sequence and Pfam parquet files, without changing their seed MSAs or HMMs. It recomputes representatives, structures, annotations, Foldseek hits, domain architectures and (with `mgnprotein_db_config`) biomes, and publishes a delta database, `db/<sample>_update.sqlite3`, holding only the successfully updated families (built by `--mode post_update_mgnifams_update_db` from the published outdir). The pipeline never modifies the production DB: `assets/merge_update_delta.sql` merges the delta in a single transaction, and the README documents it. New parameters: `--parquet_chunks`, `--hmm_chunk_size`, `--run_alphafold2`, `--colabfold_params_path`, `--af2_max_msa_seqs`, `--af2_num_recycles`.
  - `--run_alphafold2 true` also predicts each representative with ColabFold from its family full MSA (GPU). These predictions are published under `structures/alphafold2/`.
  - New `--mode post_update_mgnifams_update_db`: builds the `update_mgnifams` delta database from its outdir (`update_families/update_info.csv` records what the update computed), so the update run itself does no sqlite work.
  - New `mgnifam.seed_size` column: the number of sequences in the family seed MSA.
  - New `mgnifam.hmm_length` column: the number of HMM match states, equal to the consensus length.
  - New `mgnifam` search flags `has_pfam`, `has_funfam`, `has_model_pfam` and `has_structure` (has a Foldseek hit). `init_mgnifams_db` fills them, creates the website indexes and runs `ANALYZE` (`assets/finalize_db.sql`); `assets/merge_update_delta.sql` recomputes the flags of the updated families.
  - `mgnifam_folds` gains the Foldseek TM-scores `aln_tmscore`, `q_tmscore` and `t_tmscore`.
  - `assets/migrate_schema_seed_size_tmscores.sql` adds these columns to an existing database and backfills `seed_size` and `hmm_length`. Run it once, before the first `update_mgnifams` merge.
- `update_mgnifams` writes the new release's `update_families/families.tsv.gz` from the new samplesheet column `mgnifams_families` (the previous release's FTP `families.tsv.gz`).

### `Changed`

- `conf/base.config` sets resources for every process. The ones added (the update mode, AlphaFold2, the delta DB build, the aliased `HMMER_HMMSEARCH`, `MULTIQC`, `ARIA2`) are estimates that have not yet run at MGnify Proteins scale.
- [#67](https://github.com/EBI-Metagenomics/mgnifams/pull/67) - nf-core modules update:
  - nf-core modules updated (new versions under `Dependencies`), plus newer `hhsuite`, `hmmer/hmmsearch`, `find/concatenate` and `pigz/uncompress` containers. These modules report their versions through topic channels.
  - `update_families/family_metadata.csv` keeps the mgnifam 4.0.0 `update_families` header: `family_id,converged,seed_msa_size,full_msa_size,rep_protein,rep_region,rep_length,consensus_length,rep_sequence,consensus_sequence` (`converged` and `seed_msa_size` are empty, as `--skip_refine` builds no seed).
- [#64](https://github.com/EBI-Metagenomics/mgnifams/pull/64) - Outputs, parameters and DB layout:
  - **Breaking:** the samplesheet `sample` must match `^[A-Za-z0-9._-]+$`.
  - `hhdb_path` is required only by the default `run_mgnifams_pipeline` mode.
  - `generate_families/families/family_metadata.csv` now starts with a header row: `family_id,full_msa_size,protein,region,length,sequence,consensus,converged`.
  - `mgnifam` columns in new databases are grouped by topic (family, representative, HMM, structure, composition, flags, blobs). Databases migrated with `assets/migrate_schema_seed_size_tmscores.sql` keep the new columns last, so select columns by name, not `SELECT *` position.
  - `mgnifam_pfams` / `mgnifam_funfams` keep only domains whose i-Evalue (the exported `e_value`) is at most `--hmmsearch_evalue_cutoff`; `hmmsearch -E` filters whole sequences only, so weak domains of a significant sequence got through. The raw `domtbl` outputs are unchanged.
  - Foldseek outputs (`annotation/structures/foldseek/pdb.m8`, `all_hits.tsv`) now start with a header row (`--format-mode 4`): `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore` (tab-separated).
  - `QUERY_MGNPROTEIN_DB` outputs (`biome_mapping.tsv`, `pfam_mapping.tsv`, `query_results/`) are published under `mgnprotein_db_query/` instead of the outdir root.
  - `IMPORT_QUERIES` imports `mgnifam.csv` by header name instead of column position, with journaling and fsync off for the throwaway build file. Only the finished database is published (`INIT_SQLITE` / `IMPORT_QUERIES` outputs no longer are).
- `foldseek/easysearch` bumped to foldseek 10.941cd33. New `--foldseek_pdb_db` (default `pdb`) and `--foldseek_alphafold_db` (default `afdb50`, was hardcoded `alphafold`) name the database subdirectories under `--foldseek_db_path`.

### `Fixed`

- `-stub-run` of the default mode with `--fasta_input_mode` failed: the `EXTRACT_UNIQUE_CLUSTER_REPS` stub called `python` in an image without it, so its `versions.yml` held a shell error.
- `RUN_ESMFOLD_CPU` (the CPU fallback for long sequences and CUDA OOM failures) inherited the `process_gpu` label of `RUN_ESMFOLD` and requested a GPU with `-profile gpu`. Its `accelerator` is now cleared; SLURM site configs should request GPUs with `clusterOptions = { task.accelerator ? '--gres=gpu:1' : '' }` under `withLabel: process_gpu`.
- [#64](https://github.com/EBI-Metagenomics/mgnifams/pull/64) - Bugs found while building the update mechanism:
  - [#62](https://github.com/EBI-Metagenomics/mgnifams/issues/62) - `mgnifam_pfams` / `mgnifam_funfams` `e_value` and `score` held the full-sequence values. They now hold the per-domain i-Evalue and domain score. Rows in existing databases keep the old values until their families are updated.
  - [#58](https://github.com/EBI-Metagenomics/mgnifams/issues/58), [#59](https://github.com/EBI-Metagenomics/mgnifams/issues/59) - Domain architectures placed MGnifam domains of `<mgyp>/<start>-<end>` members at the wrong start, and ordered Pfam domains by HMM position instead of their position on the protein.
  - [#61](https://github.com/EBI-Metagenomics/mgnifams/issues/61) - `init_mgnifams_db` stored empty strings instead of NULL, and failed when an optional table CSV was missing.
  - [#60](https://github.com/EBI-Metagenomics/mgnifams/issues/60), [#63](https://github.com/EBI-Metagenomics/mgnifams/issues/63) - `-stub-run` failed in `EXTRACT_ESMFOLD_SCORES`, `PARSE_CIF` and the database modes.

### `Dependencies`

| Tool                  | Previous version | New version |
| --------------------- | ---------------- | ----------- |
| mmseqs2               | 17.b804f         | 18.8cc5c    |
| foldseek              | 9.427df8a        | 10.941cd33  |
| seqkit                | 2.9.0            | 2.13.0      |
| multiqc               | 1.29             | 1.35        |
| aria2                 | 1.36.0           | 1.37.0      |
| mgnifam               |                  | 4.0.0       |
| pyarrow               |                  | 21.0.0      |
| proteinfold_colabfold |                  | 2.0.0       |

### `Internal`

- [#64](https://github.com/EBI-Metagenomics/mgnifams/pull/64) - `prek` pre-commit hooks (`prek run --all-files`); pipeline-owned files are `nextflow lint` clean.

## v2.0.0 - [2026/04/14]

### `Added`

- [#52](https://github.com/EBI-Metagenomics/mgnifams/pull/52) - Added pfam annotations for family representatives.
- [#50](https://github.com/EBI-Metagenomics/mgnifams/pull/50) - Added support for local deeptmhmm executions.
- [#48](https://github.com/EBI-Metagenomics/mgnifams/pull/48) - Added update_db workflow.
- [#47](https://github.com/EBI-Metagenomics/mgnifams/pull/47) - Added init_db workflow.
- [#45](https://github.com/EBI-Metagenomics/mgnifams/pull/45) - Added end-to-end nf-test.
- [#41](https://github.com/EBI-Metagenomics/mgnifams/pull/41) - Added funfams annotations for family representatives.
- [#37](https://github.com/EBI-Metagenomics/mgnifams/pull/37) - Added various custom MultiQC reports; family metadata, discarded clusters, family similarities.
- Updated sqlite database subworkflow, for smooth subsequent version releases.

### `Fixed`

[#54](https://github.com/EBI-Metagenomics/mgnifams/pull/54) - ESM-1b max is 1022 residues; truncated longer sequences to avoid CRF mask errors in DeepTMHMM
[#37](https://github.com/EBI-Metagenomics/mgnifams/pull/37) - MGYP protein slices for MGnifams seed MSAs, full MSAs, reps fasta and metadata are now calculated properly in the `generate_families.py` script.

### `Changed`

- [#56](https://github.com/EBI-Metagenomics/mgnifams/pull/56) - Parallelized and optimized the `parse_domains.py` script for large-scale.
- [#55](https://github.com/EBI-Metagenomics/mgnifams/pull/55) - Parallelized and optimized the `query_mgnprotein_db.py` script for large-scale.
- [#53](https://github.com/EBI-Metagenomics/mgnifams/pull/53) - Parallelized and optimized the `pool_nonredundant_families.py` script for large-scale.
- [#42](https://github.com/EBI-Metagenomics/mgnifams/pull/42) - Replaced in house `hhsuite/reformat`, `hhsuite/hhblits` and `hhsuite/hhsearch` modules with the nf-core ones.
- [#39](https://github.com/EBI-Metagenomics/mgnifams/pull/39) - Further optimizations for MGnifams main algorithm, removing unnecessary I/O operations, by @althonos
- [#37](https://github.com/EBI-Metagenomics/mgnifams/pull/37) - Added protein set Jaccard score calculation in `identify_redundant_fams.py`, for better estimation of family similarities, in addition to `hmmsearch` among families.
  Pipeline `params` moved to root `main.nf` file, and are being passed downstream to subworkflows.
- [#31](https://github.com/EBI-Metagenomics/mgnifams/pull/31) - Swapped all subprocess calls of `generate_families.py` to cythonised lib versions (pyfamsa, pyhmmer, pytrimal).
  Benchmark results: CPU usage decrease 37.6% - Memory decrease 0.36% - Job duration decrease 37.5% - I/O read decrease 90% - I/O write decrease 77.7%

### `Dependencies`

| Tool      | Previous version | New version |
| --------- | ---------------- | ----------- |
| pytrimal  |                  | 0.8.2       |
| pyhmmer   |                  | 0.11.1      |
| pyfamsa   |                  | 0.6.0       |
| biopython |                  | 1.85        |
| pandas    |                  | 2.3.2       |
| numpy     |                  | 2.3.2       |

### `Removed`

- [#42](https://github.com/EBI-Metagenomics/mgnifams/pull/42) - Removed obsolete `reformat_msa` subworkflow, since sequence regions are now calculated inside the MGnifam algorithm.

## v1.0.0 - [2024/10/01]

Initial release of ebi-metagenomics/mgnifams.

### `Added`

- Amino acid sequence clustering (mmseqs)
- Multiple sequence alignment (mafft)
- Hidden Markov Model generation (hmmer)
- Between families redundancy removal (hhsuite)
- In-family sequence redundancy removal (esl-weight)
- Flag transmembrane families (deeptmhmm)
- Model annotation against Pfams (hhsuite)
- Structure prediction (esmfold)
- Homolog detection via structural comparison (foldseek)
