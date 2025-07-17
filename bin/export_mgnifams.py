#!/usr/bin/env python3

import argparse
import pandas as pd

def write_mgnifam_csv(metadata, structure_scores, composition, outfile):
    mgnifam_headers = ['id', 'full_size', 'protein_rep', 'rep_region', 'rep_length', 'converged', \
                        'plddt', 'ptm', 'helix_percent','strand_percent','coil_percent', \
                        'rep_sequence', 'consensus', 'seed_msa_blob', 'hmm_blob', 'rf_blob',\
                        'cif_blob', 'biome_blob', 'domain_blob', 's4pred_blob' ]

    df1 = pd.read_csv(metadata, header=None)
    df1.columns = ['id', 'full_size', 'protein_rep', 'rep_region', 'rep_length',
                    'rep_sequence', 'consensus', 'converged']

    df2 = pd.read_csv(structure_scores)

    merged = pd.merge(df1, df2, on='id', how='left')

    df3 = pd.read_csv(composition)
    merged = pd.merge(merged, df3, on='id', how='left')

    # Map columns to final header, fill missing ones
    for col in mgnifam_headers:
        if col not in merged.columns:
            merged[col] = ''

    # Reorder columns
    merged = merged[mgnifam_headers]

    merged.to_csv(outfile, index=False, header=mgnifam_headers)

def main():
    parser = argparse.ArgumentParser(description="Export MGnifams sql-ready table CSV files.")
    parser.add_argument("--metadata", help="Generated families metadata mqc CSV")
    parser.add_argument("--structure_scores", help="Tertiary prediction structure scores (plddit, ptm)")
    parser.add_argument("--composition", help="Predicted compositional features --helix, strand or coil")
    parser.add_argument("--outfile", help="CSV for mgnifam table")

    args = parser.parse_args()

    write_mgnifam_csv(args.metadata, args.structure_scores, args.composition, args.outfile)

if __name__ == "__main__":
    main()
