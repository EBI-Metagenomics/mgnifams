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

### Publishing an update release

How 1.1 was published (2026-11-02), for the next `X.Y` after `X.(Y-1)`. On Codon, `/nfs/ftp` is mounted only on the
datamover nodes, so every step that touches it runs there (`sbatch -p datamover`). The FTP is a separate filesystem
from `/nfs/production`, so files are copied to it, not hard-linked.

1. **Back up** the previous release's `full_msa.tar*` outside the FTP, and check it against that release's `md5sums.txt`.
2. **Stage** the new files outside the public tree:
   `sbatch -p datamover build_update_release.sh <update outdir> <FTP>/releases/X.(Y-1) <stage dir>`. It builds the
   complete `full_msa.tar` (the previous alignments of the `not_updated` families), copies `families.tsv.gz` and
   `pipeline_info/`, and writes `md5sums.txt`. It stops on any count mismatch.
3. **Fill** `releases/X.Y/release_manifest.json`, `README.md` and the `RELEASES.tsv` row in this repository, from the
   staging log (md5s and sizes) and the run's `pipeline_info/`.
4. **Copy into a hidden folder** `<FTP>/releases/.X.Y_tmp`: hard-link `mgnifams_hmm.lib.gz` and `seed_msa.tar.gz` from
   `releases/X.(Y-1)/`, copy the staged files and the metadata from this repository, then run `md5sum -c md5sums.txt`
   there.
5. **Go live:** rename `.X.Y_tmp` to `X.Y`, swap the link in one step
   (`ln -s releases/X.Y .cr && mv -T .cr current_release`), and replace `RELEASES.tsv` and the top-level `README.md`
   (without this repository's "Draft" note).
6. **Retire** the previous `full_msa.tar*`: remove it and its line in that release's `md5sums.txt` (the checksum stays
   in its `release_manifest.json`). Then delete the stage dir.

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
| `rep_region`               | representative region, `start-end` (1-based, inclusive); `-` = the whole protein (1.0 only)   |

MGnify Proteins releases are cited by their numbered folder (`peptide_database/2026_07/`), never by
`current_release`.
