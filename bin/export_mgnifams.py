#!/usr/bin/env python3

import argparse
import pandas as pd
import os


def write_mgnifam_csv(metadata, structure_scores, composition, tm_composition, outfile, seed_sizes=""):
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

    # family_metadata.csv header (family_id,full_msa_size,protein,...) renamed to the mgnifam column names
    df1 = pd.read_csv(metadata)
    df1.columns = [
        "id",
        "full_size",
        "protein_rep",
        "rep_region",
        "rep_length",
        "rep_sequence",
        "consensus",
        "converged",
    ]

    # HMM length = one consensus residue per match state
    df1["hmm_length"] = df1["consensus"].astype("string").str.len().astype("Int64")

    df2 = pd.read_csv(structure_scores)
    merged = pd.merge(df1, df2, on="id", how="left")

    df3 = pd.read_csv(composition)
    merged = pd.merge(merged, df3, on="id", how="left")

    # Optional: TM composition
    if tm_composition and os.path.isfile(tm_composition) and os.path.getsize(tm_composition) > 0:
        df4 = pd.read_csv(tm_composition)
        merged = pd.merge(merged, df4, on="id", how="left")
    else:
        print("Note: TM composition file not provided or empty. Skipping...")

    # Optional: seed MSA sizes (id,seed_size); absent in update mode, where seeds do not change
    if seed_sizes and os.path.isfile(seed_sizes) and os.path.getsize(seed_sizes) > 0:
        merged = pd.merge(merged, pd.read_csv(seed_sizes), on="id", how="left")
        merged["seed_size"] = merged["seed_size"].astype("Int64")

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
    parser.add_argument("--structure_scores", required=True, help="Tertiary prediction structure scores (plddt, ptm)")
    parser.add_argument("--composition", required=True, help="Predicted compositional features --helix, strand or coil")
    parser.add_argument(
        "--tm_composition", default="", help="Predicted transmembrane features --inside, membrane or outside"
    )
    parser.add_argument("--seed_sizes", default="", help="Optional CSV with id,seed_size (seed MSA sequence counts)")
    parser.add_argument("--outfile", required=True, help="CSV for mgnifam table")

    args = parser.parse_args()

    write_mgnifam_csv(
        args.metadata, args.structure_scores, args.composition, args.tm_composition, args.outfile, args.seed_sizes
    )


if __name__ == "__main__":
    main()
