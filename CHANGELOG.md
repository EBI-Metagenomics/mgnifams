# ebi-metagenomics/mgnifams: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [yyyy/mm/dd]

### `Added`

- Update sqlite database subworkflow, for smooth subsequent version releases.

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
