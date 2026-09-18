# ebi-metagenomics/mgnifams: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.0dev - [unreleased]

### `Added`

- New `--mode update_mgnifams`: refreshes existing families against new MGnify proteins, read from the MGnify proteins sequence and Pfam parquet files, without changing their seed MSAs or HMMs. It recomputes representatives, structures, annotations, Foldseek hits, domain architectures and (with `mgnprotein_db_config`) biomes, and publishes a delta database, `db/<sample>_update.sqlite3`, holding only the successfully updated families. The pipeline never modifies the production DB: `assets/merge_update_delta.sql` merges the delta in a single transaction, and the README documents it. New parameters: `--parquet_chunks`, `--hmm_chunk_size`, `--run_alphafold2`, `--colabfold_params_path`, `--af2_max_msa_seqs`, `--af2_num_recycles`.
  - `--run_alphafold2 true` also predicts each representative with ColabFold from its family full MSA (GPU). These predictions are published under `structures/alphafold2/`.
- New `mgnifam.seed_size` column: the number of sequences in the family seed MSA.
- `mgnifam_folds` gains the Foldseek TM-scores `aln_tmscore`, `q_tmscore` and `t_tmscore`.
- `assets/migrate_schema_seed_size_tmscores.sql` adds both of these to an existing database and backfills `seed_size`. Run it once, before the first `update_mgnifams` merge.

### `Changed`

- **Breaking:** the samplesheet `sample` must match `^[A-Za-z0-9._-]+$`.
- `hhdb_path` is required only by the default `run_mgnifams_pipeline` mode.

### `Fixed`

- [#62](https://github.com/EBI-Metagenomics/mgnifams/issues/62) - `mgnifam_pfams` / `mgnifam_funfams` `e_value` and `score` held the full-sequence values. They now hold the per-domain i-Evalue and domain score. Rows in existing databases keep the old values until their families are updated.
- [#58](https://github.com/EBI-Metagenomics/mgnifams/issues/58), [#59](https://github.com/EBI-Metagenomics/mgnifams/issues/59) - Domain architectures placed MGnifam domains of `<mgyp>/<start>-<end>` members at the wrong start, and ordered Pfam domains by HMM position instead of their position on the protein.
- [#61](https://github.com/EBI-Metagenomics/mgnifams/issues/61) - `init_mgnifams_db` stored empty strings instead of NULL, and failed when an optional table CSV was missing.
- [#60](https://github.com/EBI-Metagenomics/mgnifams/issues/60), [#63](https://github.com/EBI-Metagenomics/mgnifams/issues/63) - `-stub-run` failed in `EXTRACT_ESMFOLD_SCORES`, `PARSE_CIF` and the database modes.

### `Internal`

- `prek` pre-commit hooks (`prek run --all-files`); pipeline-owned files are `nextflow lint` clean.

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
