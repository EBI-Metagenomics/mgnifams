# ebi-metagenomics/mgnifams: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [yyyy/mm/dd]

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

[#37](https://github.com/EBI-Metagenomics/mgnifams/pull/37) - MGYP protein slices for MGnifams seed MSAs, full MSAs, reps fasta and metadata are now calculated properly in the `generate_families.py` script.

### `Changed`

- [#53](https://github.com/EBI-Metagenomics/mgnifams/pull/53) - Parallelized and optimized the `pool_nonredundant_families.py` script for large-scale.
- [#42](https://github.com/EBI-Metagenomics/mgnifams/pull/42) - Replaced in house `hhsuite/reformat`, `hhsuite/hhblits` and `hhsuite/hhsearch` modules with the nf-core ones.
- [#39](https://github.com/EBI-Metagenomics/mgnifams/pull/39) - Further optimizations for MGnifams main algorithm, removing unnecessary I/O operations, by @althonos
- [#37](https://github.com/EBI-Metagenomics/mgnifams/pull/37) - Added protein set Jaccard score calculation in `identify_redundant_fams.py`, for better estimation of family similarities, in addition to `hmmsearch` among families.
  Pipeline `params` moved to root `main.nf` file, and are being passed downstream to subworkflows.
- [#31](https://github.com/EBI-Metagenomics/mgnifams/pull/31) - Swapped all subprocess calls of `generate_families.py` to cythonised lib versions (pyfamsa, pyhmmer, pytrimal).
Benchmark results: CPU usage decrease 37.6% - Memory decrease 0.36% - Job duration decrease 37.5% - I/O read decrease 90% - I/O write decrease 77.7%

### `Dependencies`

| Tool       | Previous version | New version |
| ---------- | ---------------- | ----------- |
| pytrimal   |                  | 0.8.2       |
| pyhmmer    |                  | 0.11.1      |
| pyfamsa    |                  | 0.6.0       |
| biopython  |                  | 1.85        |
| pandas     |                  | 2.3.2       |
| numpy      |                  | 2.3.2       |

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
