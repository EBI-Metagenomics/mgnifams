# MGnifams FTP

> Draft of the layout of <https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnifams/>. Only the READMEs and the small
> metadata files are kept in this repository; the archives are listed but not committed.

MGnifams are protein families built from metagenomics-derived protein sequences of
[MGnify Proteins](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/) by the
[MGnifams pipeline](https://github.com/EBI-Metagenomics/mgnifams).

## Releases

MGnifams releases have their own numbers, `MAJOR.MINOR`, separate from the MGnify Proteins releases (`YYYY_MM`).
Each release records the MGnify Proteins release it was built on.

| Change                                                           | Bump    | HMM library           |
| ---------------------------------------------------------------- | ------- | --------------------- |
| Family members refreshed against a newer MGnify Proteins release | `MINOR` | unchanged (same file) |
| New families added, or seeds / HMMs rebuilt                      | `MAJOR` | changed               |

The HMMs within one `MAJOR` release are identical, so `hmmsearch` / `hmmscan` results against any `1.x` library stay valid.

[`RELEASES.tsv`](RELEASES.tsv) lists every release, with the same columns as the `release_history` table of the MGnifams database:

| Release | Type    | MGnify Proteins | Families |
| ------- | ------- | --------------- | -------- |
| 1.0     | initial | 2024_04         | 35,459   |
| 1.1     | update  | 2026_07         | 35,459   |

## Layout

```
mgnifams/
|-- README                    this file
|-- RELEASES.tsv              one row per release
|-- current_release -> releases/1.1
`-- releases/
    |-- 1.0/
    |   |-- README
    |   |-- release_manifest.json
    |   |-- md5sums.txt
    |   |-- families.tsv.gz
    |   |-- mgnifams_hmm.lib.gz
    |   |-- seed_msa.tar.gz
    |   `-- pipeline_info/        final parameters (without local paths) and software versions of the run
    `-- 1.1/
        |-- README                includes the changes since 1.0
        |-- release_manifest.json
        |-- md5sums.txt
        |-- families.tsv.gz
        |-- mgnifams_hmm.lib.gz   same file as 1.0 (hard link, same md5)
        |-- seed_msa.tar.gz       same file as 1.0 (hard link, same md5)
        |-- full_msa.tar          latest release only
        `-- pipeline_info/
```

`current_release` always points to the latest release. To cite a release, or to download the same files again, use
its numbered folder, since `current_release` moves.

The files that used to sit at the top level (`mgnifams_hmm.lib.gz`, `seed_msa.tar.gz`) link to `releases/1.0/`, so
existing links still return the same files. The old top-level `full_msa.tar.gz` link is removed together with
`releases/1.0/full_msa.tar.gz` when 1.1 is published.

### Moving to this layout

From the flat layout (1.0 files at the top level):

1. Create `releases/1.0/` and move `mgnifams_hmm.lib.gz`, `seed_msa.tar.gz` and `full_msa.tar.gz` into it. Add its
   `README`, `release_manifest.json`, `families.tsv.gz`, `md5sums.txt` and `pipeline_info/`, then run `md5sum -c md5sums.txt`
   there.
2. Link `current_release -> releases/1.0`, and the three old top-level names to their files in `releases/1.0/`.
3. Replace the top-level `README` with this file, and add `RELEASES.tsv` with the 1.0 row only.

When 1.1 is published:

1. Create `releases/1.1/`: hard-link `mgnifams_hmm.lib.gz` and `seed_msa.tar.gz` from `releases/1.0/`, and add
   `full_msa.tar` and the metadata files (`families.tsv.gz` is the update run's `update_families/families.tsv.gz`).
   Check it with `md5sum -c md5sums.txt`.
2. Point `current_release` to `releases/1.1` and add the 1.1 row to `RELEASES.tsv`.
3. Remove `releases/1.0/full_msa.tar.gz` and the top-level `full_msa.tar.gz` link.

## Files

| File                    | Kept in             | Content                                                                                                         |
| ----------------------- | ------------------- | --------------------------------------------------------------------------------------------------------------- |
| `mgnifams_hmm.lib.gz`   | every release       | HMMER3 profile library, one HMM per family, named by the integer family id                                      |
| `seed_msa.tar.gz`       | every release       | per-family seed alignments, `<family_id>.fas.gz` (aligned FASTA, gaps as `-`), from which the HMMs are built    |
| `full_msa.tar`          | latest release only | per-family full alignments, `<family_id>.sto.gz` (Stockholm): every MGnify sequence recruited by the family HMM |
| `families.tsv.gz`       | every release       | one row per family with its provenance (below)                                                                  |
| `release_manifest.json` | every release       | release metadata, inputs, tool and database versions, file checksums                                            |
| `md5sums.txt`           | every release       | `md5sum -c` compatible checksums of the release files                                                           |

`full_msa.tar` is a plain tar: each member is already gzipped, so compressing the archive again would not make it
smaller. It is always complete. Families that a release did not update keep the alignment of the release that last
updated them, as recorded in `families.tsv.gz`. When a new release is published, the full MSA archive of the previous
one is removed; its checksum stays in that release's manifest.

### `families.tsv.gz`

Tab-separated, with a header row. Families are never removed or renumbered, and family ids are never reused.

| Column                     | Meaning                                                                                       |
| -------------------------- | --------------------------------------------------------------------------------------------- |
| `family_id`                | integer MGnifams family id (the HMM `NAME`)                                                   |
| `status`                   | `new` (first appears in this release), `updated` or `not_updated`                             |
| `not_updated_reason`       | why an update left the family as it was (e.g. `no hits in the new database`); empty otherwise |
| `first_release`            | release in which the family first appeared                                                    |
| `model_release`            | release that built the family's seed MSA and HMM                                              |
| `model_proteins_release`   | MGnify Proteins release of the seed MSA                                                       |
| `members_release`          | release that last computed the full MSA, representative and annotations                       |
| `members_proteins_release` | MGnify Proteins release of the full MSA members                                               |
| `full_size`                | number of sequences in the full MSA                                                           |
| `seed_size`                | number of sequences in the seed MSA                                                           |
| `rep_id`                   | representative MGnify protein (`MGYP...`)                                                     |
| `rep_region`               | representative region, `start-end` (1-based, inclusive)                                       |

MGnify Proteins releases are cited by their numbered folder (`peptide_database/2026_07/`), never by
`current_release`.
