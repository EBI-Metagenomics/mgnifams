import sys
import pandas as pd
import csv

def read_hh_hits(hh_hits_file):
    column_names = ["Fam", "Hit", "Prob", "E-value", "P-value", "Score", "SS", "Cols", "Query HMM", "Template HMM", ""]
    hh_hits      = pd.read_csv(hh_hits_file, delimiter='\t', header=None, names=column_names, skiprows=1)

    return hh_hits

def read_rep_to_fam_dict(fam_rep_mapping_file):
    rep_to_fam_dict = {}
    with open(fam_rep_mapping_file, mode='r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            key = row[1]
            value = row[0]
            rep_to_fam_dict[key] = value

    return rep_to_fam_dict

def map_and_remove_self(hh_hits, rep_to_fam_dict):
    hh_hits['Hit'] = hh_hits['Hit'].map(rep_to_fam_dict) # Do the mapping on the df
    hh_hits = hh_hits[hh_hits['Fam'] != hh_hits['Hit']]
    return hh_hits

def remove_redundant(hh_hits, fam):
    hh_hits = hh_hits[(hh_hits['Fam'] != fam) & (hh_hits['Hit'] != fam)]
    return hh_hits

def remove_pair(hh_hits, fam, hit):
    hh_hits = hh_hits[(hh_hits['Fam'] != fam) & (hh_hits['Hit'] != hit)]
    hh_hits = hh_hits[(hh_hits['Fam'] != hit) & (hh_hits['Hit'] != fam)]
    return hh_hits

def check_similarity_remove_if_redundant(hh_hits, row, \
    fams_to_export, similarity_edgelist_file, \
    redundant_threshold, similarity_threshold):
    fam       = row['Fam']
    hit       = row['Hit']
    avg_prob  = hh_hits[(hh_hits['Fam'] == fam) & (hh_hits['Hit'] == hit)]['Prob'].mean()
    
    if (avg_prob >= redundant_threshold):
        # keeping fam with better Score (e.g. sum all fam1 hit4), if same, more Cols(e.g. sum all fam1 hit4)
        avg_score1 = hh_hits[(hh_hits['Fam'] == fam) & (hh_hits['Hit'] == hit)]['Score'].mean()
        avg_score2 = hh_hits[(hh_hits['Fam'] == hit) & (hh_hits['Hit'] == fam)]['Score'].mean()
        if (avg_score1 != avg_score2): # check Score
            if (avg_score1 > avg_score2):
                fams_to_export.remove(hit)
                hh_hits = remove_redundant(hh_hits, hit)
            else:
                fams_to_export.remove(fam)
                hh_hits = remove_redundant(hh_hits, fam)
        else: # check Cols
            avg_cols1  = hh_hits[(hh_hits['Fam'] == fam) & (hh_hits['Hit'] == hit)]['Cols'].mean()
            avg_cols2  = hh_hits[(hh_hits['Fam'] == hit) & (hh_hits['Hit'] == fam)]['Cols'].mean()
            if (avg_cols1 >= avg_cols2):
                fams_to_export.remove(hit)
                hh_hits = remove_redundant(hh_hits, hit)
            else:
                fams_to_export.remove(fam)
                hh_hits = remove_redundant(hh_hits, fam)
    elif (avg_prob >= similarity_threshold):
        with open(similarity_edgelist_file, 'a') as f: 
            f.write(f"{fam},{hit},{avg_prob}\n")
        hh_hits = remove_pair(hh_hits, fam, hit)

    return hh_hits, fams_to_export

def write_non_redundant_fam_ids(fams_to_export):
    with open(non_redundant_fam_ids_file, 'w') as f:
        for id in fams_to_export:
            f.write(f"{id}\n")

def export_non_redundant_family_ids(hh_hits_file, fam_rep_mapping_file, \
    non_redundant_fam_ids_file, similarity_edgelist_file, \
    redundant_threshold=100, similarity_threshold=50): # TODO redundant_threshold=95, testing

    hh_hits         = read_hh_hits(hh_hits_file)
    fams_to_export  = hh_hits['Fam'].unique().tolist() # This must be done here, before removing self-hits (some fams might have only self-hits)
    rep_to_fam_dict = read_rep_to_fam_dict(fam_rep_mapping_file)
    hh_hits         = map_and_remove_self(hh_hits, rep_to_fam_dict)
    
    i = 0
    while len(hh_hits) > i:
        row = hh_hits.iloc[i]
        hh_hits, fams_to_export = check_similarity_remove_if_redundant(hh_hits, row, \
            fams_to_export, similarity_edgelist_file, redundant_threshold, similarity_threshold)
        i += 1

    write_non_redundant_fam_ids(fams_to_export)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_redundant.py <hh_hits> <fam_rep_mapping> <output_fam_ids>")
        sys.exit(1)

    non_redundant_fam_ids_file = "non_redundant_fam_ids.txt"
    similarity_edgelist_file   = "similarity_edgelist.csv"
    export_non_redundant_family_ids(sys.argv[1], sys.argv[2], non_redundant_fam_ids_file, similarity_edgelist_file)
