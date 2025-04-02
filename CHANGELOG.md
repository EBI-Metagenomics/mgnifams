# ebi-metagenomics/mgnifams: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [yyyy/mm/dd]

### `Added`

- Update sqlite database subworkflow, for smooth subsequent version releases.

### `Changed`

- [#31](https://github.com/EBI-Metagenomics/mgnifams/pull/31) - Swapped all subprocess calls of `generate_families.py` to cythonised lib versions (pyfastx, pyfamsa, pyhmmer, pytrimal).
Benchmark results: CPU usage decrease 37.6% - Memory decrease 0.36% - Job duration decrease 37.5% - I/O read decrease 90% - I/O write decrease 77.7%

### `Dependencies`

| Tool     | Previous version | New version |
| -------- | ---------------- | ----------- |
| pytrimal |                  | 0.7.0       |
| pyhmmer  |                  | 0.11.0      |
| pyfamsa  |                  | 0.5.3.post1 |
| pyfastx  |                  | 2.2.0       |

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
