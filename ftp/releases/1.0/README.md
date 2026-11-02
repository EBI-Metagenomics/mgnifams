# MGnifams 1.0

Initial release: 35,459 families built from [MGnify Proteins 2024_04](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2024_04/).

| File                    | Content                                                                    |
| ----------------------- | -------------------------------------------------------------------------- |
| `mgnifams_hmm.lib.gz`   | HMMER3 profile library, one HMM per family, named by the integer family id |
| `seed_msa.tar.gz`       | per-family seed alignments, `<family_id>.fas.gz`                           |
| `families.tsv.gz`       | one row per family with its provenance (all `new` in this release)         |
| `release_manifest.json` | inputs, tool and database versions, file checksums                         |
| `md5sums.txt`           | checksums of the files above                                               |
| `pipeline_info/`        | final parameters (without local paths) and software versions of the run    |

`full_msa.tar.gz` was removed when 1.1 was published: only the latest release keeps its full MSAs. Its checksum
stays in `release_manifest.json`. For the current full MSAs, see [`current_release`](../../current_release).

The HMMs and seed MSAs are the same in every `1.x` release.

## Known issues

These come from the pipeline version that built 1.0 (`v2.0.0dev`); later versions fixed them. `families.tsv.gz`
reports the published data as it is.

| Families | Issue                                               |
| -------- | --------------------------------------------------- |
| 544      | seed MSA with a single sequence                     |
| 3,126    | seed MSA above the 2,000-sequence cap (up to 2,086) |
| 2        | full MSA smaller than the seed MSA (12187, 12923)   |
| 397      | no representative region (`rep_region` is `-`)      |

1.1 recomputes the representatives of the families it updates. The seed MSAs stay as they are in every `1.x` release.

Format details are in the [top-level README](../../README.md).
