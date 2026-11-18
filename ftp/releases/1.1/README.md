# MGnifams 1.1

Update of the 1.0 families against [MGnify Proteins 2026_07](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/),
produced by the `update_mgnifams` mode of the [MGnifams pipeline](https://github.com/EBI-Metagenomics/mgnifams).

| File                    | Content                                                                   |
| ----------------------- | ------------------------------------------------------------------------- |
| `mgnifams_hmm.lib.gz`   | HMMER3 profile library, **the same file as 1.0**                          |
| `seed_msa.tar.gz`       | per-family seed alignments, **the same file as 1.0**                      |
| `full_msa.tar`          | per-family full alignments, `<family_id>.sto.gz`, for all 35,459 families |
| `families.tsv.gz`       | one row per family with its provenance and update status                  |
| `clans.tsv.gz`          | clan membership: 4,363 clans of 20,964 families (added 2026-11-09)        |
| `clan_edges.tsv.gz`     | residue overlap of similar family pairs, from which the clans are built   |
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
- **New:** clans, for all 35,459 families (below). 1.0 has no clan file.
- **Kept at their 1.0 version:** families the update did not refresh (`status` = `not_updated`, with the reason in
  `not_updated_reason`). Their alignment in `full_msa.tar` is the 1.0 one (`members_proteins_release` = `2024_04`).

| Status        | Families |
| ------------- | -------- |
| `updated`     | 34,419   |
| `not_updated` | 1,040    |

A family is `not_updated` when its 1.0 HMM, searched against 2026_07, did not give a usable new full MSA:

| `not_updated_reason`                         | Families | Meaning in an update                                                           |
| -------------------------------------------- | -------- | ------------------------------------------------------------------------------ |
| `family representative length too small`     | 788      | the new representative (row 0 of the new full MSA) is shorter than 75 residues |
| `low complexity model - confounding cluster` | 215      | the HMM found hits in 2026_07, but none covered at least 90% of the model      |
| `family representative length too large`     | 37       | the new representative is longer than 2,000 residues                           |

The name `low complexity model - confounding cluster` comes from family generation, where it has the same test (no
hit covering 90% of the model); in an update it does not mean that the model is low complexity.

## Clans

Computed with `mgnifam compute_clans` 4.1.0 ([mgnifam](https://github.com/vagkaratzas/mgnifam)) from the full MSA
members of all 35,459 families: the 1.1 members of the 34,419 updated families and the 1.0 members of the 1,040 not
updated. For two families sharing a protein, the residue Jaccard is the number of (protein, position) pairs covered by
both families' member regions over the number covered by either. A clan is a connected component of the graph whose
edges are the family pairs with residue Jaccard >= 0.5; 14,495 families are in no clan and are not listed. Clan ids
(`1.1_clan_<n>`, numbered by size, largest first) are not stable across releases. Columns are in the
[top-level README](../../README.md).

## Known limitations

- **Representative choice.** The representative is the best-scoring domain of the top-scoring sequence. When that
  sequence is a long repeat protein, its best domain can be a fragment of the model, or one long region spanning many
  repeats: 19.1% of the updated representatives are shorter than 90% of the model length (17.6% in 1.0), and most of
  the length-based `not_updated` families above come from this. It will change in the next family remake
  ([mgnifam#23](https://github.com/vagkaratzas/mgnifam/issues/23)).
- **Full MSAs include partial matches.** As in 1.0, the full MSA keeps every domain with E-value <= 0.001, whatever
  its length; the 90% length requirement only applies while the family is built.
- **`rep_region` `-`.** `rep_region` = `-` means the representative is the whole protein, in updated and
  not-updated families alike.

Format details are in the [top-level README](../../README.md).
