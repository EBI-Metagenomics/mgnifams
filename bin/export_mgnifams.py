#!/usr/bin/env python3

import argparse
import pandas as pd
import os


def write_mgnifam_csv(metadata, structure_scores, composition, tm_composition, outfile):
    mgnifam_headers = [
        "id",
        "full_size",
        "seed_size",
        "converged",
        "protein_rep",
        "rep_region",
        "rep_length",
        "rep_sequence",
        "consensus",
        "hmm_length",
        "plddt",
        "ptm",
        "helix_percent",
        "strand_percent",
        "coil_percent",
        "inside_percent",
        "membrane_alpha_percent",
        "outside_percent",
        "signal_percent",
        "membrane_beta_percent",
        "periplasm_percent",
        "seed_msa_blob",
        "hmm_blob",
        "rf_blob",
        "cif_blob",
        "biome_blob",
        "domain_blob",
        "s4pred_blob",
        "tm_blob",
    ]

    # family_metadata.csv (mgnifam generate_families / update_families: family_id,converged,seed_msa_size,
    # full_msa_size,rep_protein,rep_region,rep_length,consensus_length,rep_sequence,consensus_sequence)
    # columns renamed to the mgnifam table's column names
    df1 = pd.read_csv(metadata).rename(
        columns={
            "family_id": "id",
            "full_msa_size": "full_size",
            "seed_msa_size": "seed_size",
            "rep_protein": "protein_rep",
            "consensus_sequence": "consensus",
            "consensus_length": "hmm_length",
        }
    )

    # Absent (left empty) for ids-only metadata (annotate_not_updated)
    if "hmm_length" in df1.columns:
        df1["hmm_length"] = df1["hmm_length"].astype("Int64")
    # Empty under update_families --skip_refine (no seed built); the merge keeps prod's value
    if "seed_size" in df1.columns:
        df1["seed_size"] = df1["seed_size"].astype("Int64")

    # Optional: structure scores; absent in annotate_not_updated, where structures do not change
    if structure_scores and os.path.isfile(structure_scores) and os.path.getsize(structure_scores) > 0:
        merged = pd.merge(df1, pd.read_csv(structure_scores), on="id", how="left")
    else:
        merged = df1

    df3 = pd.read_csv(composition)
    merged = pd.merge(merged, df3, on="id", how="left")

    # Optional: TM composition
    if tm_composition and os.path.isfile(tm_composition) and os.path.getsize(tm_composition) > 0:
        df4 = pd.read_csv(tm_composition)
        merged = pd.merge(merged, df4, on="id", how="left")
    else:
        print("Note: TM composition file not provided or empty. Skipping...")

    # Ensure all expected columns are present
    for col in mgnifam_headers:
        if col not in merged.columns:
            merged[col] = ""

    # Reorder and export
    merged = merged[mgnifam_headers]
    merged.to_csv(outfile, index=False, header=mgnifam_headers)


def main():
    parser = argparse.ArgumentParser(description="Export MGnifams sql-ready table CSV files.")
    parser.add_argument("--metadata", required=True, help="family_metadata.csv (with header)")
    parser.add_argument(
        "--structure_scores", default="", help="Optional tertiary prediction structure scores (plddt, ptm)"
    )
    parser.add_argument("--composition", required=True, help="Predicted compositional features --helix, strand or coil")
    parser.add_argument(
        "--tm_composition", default="", help="Predicted transmembrane features --inside, membrane or outside"
    )
    parser.add_argument("--outfile", required=True, help="CSV for mgnifam table")

    args = parser.parse_args()

    write_mgnifam_csv(args.metadata, args.structure_scores, args.composition, args.tm_composition, args.outfile)


if __name__ == "__main__":
    main()
