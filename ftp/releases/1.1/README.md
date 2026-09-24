# MGnifams 1.1

Update of the 1.0 families against [MGnify Proteins 2026_07](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/),
produced by the `update_mgnifams` mode of the [MGnifams pipeline](https://github.com/EBI-Metagenomics/mgnifams).

| File                    | Content                                                                   |
| ----------------------- | ------------------------------------------------------------------------- |
| `mgnifams_hmm.lib.gz`   | HMMER3 profile library, **the same file as 1.0**                          |
| `seed_msa.tar.gz`       | per-family seed alignments, **the same file as 1.0**                      |
| `full_msa.tar`          | per-family full alignments, `<family_id>.sto.gz`, for all 35,459 families |
| `families.tsv.gz`       | one row per family with its provenance and update status                  |
| `release_manifest.json` | inputs, tool and database versions, file checksums                        |
| `md5sums.txt`           | checksums of the files above                                              |
| `pipeline_info/`        | final parameters (without local paths) and software versions of the run   |

## Changes since 1.0

- **Unchanged:** the family set, the family ids, the seed MSAs and the HMMs.
- **Recomputed** for every family the update succeeded on (`status` = `updated`), from its 2026_07 members:
  - the full MSA and family size;
  - the representative sequence and its structure;
  - the annotations of the representative: Pfam, FunFams, secondary structure and Foldseek hits;
  - the domain architectures and the biomes.
- **Kept at their 1.0 version:** families the update did not refresh (`status` = `not_updated`, with the reason in
  `not_updated_reason`). Their alignment in `full_msa.tar` is the 1.0 one (`members_proteins_release` = `2024_04`).

| Status        | Families |
| ------------- | -------- |
| `updated`     | TBD      |
| `not_updated` | TBD      |

Format details are in the [top-level README](../../README.md).
