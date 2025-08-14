#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def write_mgnifam_csv(metadata, structure_scores, composition, tm_composition, outfile):
    mgnifam_headers = [
        'id', 'full_size', 'protein_rep', 'rep_region', 'rep_length', 'converged',
        'plddt', 'ptm', 'helix_percent', 'strand_percent', 'coil_percent',
        'inside_percent', 'membrane_alpha_percent', 'outside_percent',
        'signal_percent', 'membrane_beta_percent', 'periplasm_percent',
        'rep_sequence', 'consensus', 'seed_msa_blob', 'hmm_blob', 'rf_blob',
        'cif_blob', 'biome_blob', 'domain_blob', 's4pred_blob', 'tm_blob'
    ]

    df1 = pd.read_csv(metadata, header=None)
    df1.columns = ['id', 'full_size', 'protein_rep', 'rep_region', 'rep_length',
                   'rep_sequence', 'consensus', 'converged']

    df2 = pd.read_csv(structure_scores)
    merged = pd.merge(df1, df2, on='id', how='left')

    df3 = pd.read_csv(composition)
    merged = pd.merge(merged, df3, on='id', how='left')

    # Optional: TM composition
    if tm_composition and os.path.isfile(tm_composition) and os.path.getsize(tm_composition) > 0:
        df4 = pd.read_csv(tm_composition)
        merged = pd.merge(merged, df4, on='id', how='left')
    else:
        print("Note: TM composition file not provided or empty. Skipping...")

    # Ensure all expected columns are present
    for col in mgnifam_headers:
        if col not in merged.columns:
            merged[col] = ''

    # Reorder and export
    merged = merged[mgnifam_headers]
    merged.to_csv(outfile, index=False, header=mgnifam_headers)

def main():
    parser = argparse.ArgumentParser(description="Export MGnifams sql-ready table CSV files.")
    parser.add_argument("--metadata", required=True, help="Generated families metadata mqc CSV")
    parser.add_argument("--structure_scores", required=True, help="Tertiary prediction structure scores (plddt, ptm)")
    parser.add_argument("--composition", required=True, help="Predicted compositional features --helix, strand or coil")
    parser.add_argument("--tm_composition", default="", help="Predicted transmembrane features --inside, membrane or outside")
    parser.add_argument("--outfile", required=True, help="CSV for mgnifam table")

    args = parser.parse_args()

    write_mgnifam_csv(args.metadata, args.structure_scores, args.composition, args.tm_composition, args.outfile)

if __name__ == "__main__":
    main()
