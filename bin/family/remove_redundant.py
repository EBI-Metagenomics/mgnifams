import sys
import pandas as pd

def export_non_redundant_family_ids(hh_hits_file, fam_rep_mapping_file, non_redundant_fam_ids_file, \
    redundant_threshold=100, similarity_threshold=50): # TODO redundant_threshold=95, testing

    column_names    = ["Fam", "Hit", "Prob", "E-value", "P-value", "Score", "SS", "Cols", "Query HMM", "Template HMM", ""]
    hh_hits         = pd.read_csv(hh_hits_file, delimiter='\t', header=None, names=column_names)
    print(hh_hits)
    fam_rep_mapping = pd.read_csv(fam_rep_mapping_file, header=None)
    print(fam_rep_mapping)
    # TODO, keep fam with better Score, if same, more Cols
    pass # TODO

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python remove_redundant.py <hh_hits> <fam_rep_mapping> <output_fam_ids>")
        sys.exit(1)

    export_non_redundant_family_ids(sys.argv[1], sys.argv[2], sys.argv[3])
