#!/usr/bin/env python3

import sys
import pandas as pd
import csv

# Initiation Functions Start ###
def initialize_outfiles(file_paths):
    """Initialize multiple files to be empty."""
    for file_path in file_paths:
        with open(file_path, "w") as file:
            pass

def create_fasta_to_length_dict(fasta_file):
    with open(log_file, 'a') as file:
        file.write("Reading family rep lengths from fasta...")

    length_dict = {}
    with open(fasta_file, 'r') as file:
        current_id = None
        current_length = 0

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    length_dict[current_id] = current_length
                current_id = line[1:]  # Get the ID (remove '>')
                current_length = 0  # Reset length for the new sequence
            else:
                current_length += len(line)

        # Adding last sequence
        if current_id is not None:
            length_dict[current_id] = current_length

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return length_dict

def read_rep_to_fam_dicts(fam_rep_mapping_file):
    with open(log_file, 'a') as file:
        file.write("Reading rep_to_fam and fam_to_rep dictionaries...")

    rep_to_fam_dict = {}
    fam_to_rep_dict = {}
    
    with open(fam_rep_mapping_file, mode='r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            rep = row[1]
            fam = row[0]
            rep_to_fam_dict[rep] = fam
            fam_to_rep_dict[fam] = rep

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return rep_to_fam_dict, fam_to_rep_dict

def read_hh_hits(hh_hits_file):
    with open(log_file, 'a') as file:
        file.write("Reading hh_hits files...")

    column_names = ["Fam", "Hit", "Prob", "E-value", "P-value", "Score", "SS", "Cols", "Query HMM", "Template HMM", ""]
    hh_hits      = pd.read_csv(hh_hits_file, delimiter='\t', header=None, names=column_names, skiprows=1)

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return hh_hits

def map_and_remove_self(hh_hits):
    with open(log_file, 'a') as file:
        file.write("Mapping reps to fams and removing self hits...")

    hh_hits['Hit'] = hh_hits['Hit'].map(rep_to_fam_dict) # Do the mapping on the df
    hh_hits = hh_hits[hh_hits['Fam'] != hh_hits['Hit']]

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return hh_hits

def keep_unique_pairs(hh_hits):
    with open(log_file, 'a') as file:
        file.write("Keeping unique pairs...")

    hh_hits = hh_hits[['Fam', 'Hit']]
    # Remove exact duplicates
    hh_hits = hh_hits.drop_duplicates()
    # Remove reverse pairs
    pairs = set()
    unique_rows = []

    for _, row in hh_hits.iterrows():
        fam_hit = (row['Fam'], row['Hit'])
        hit_fam = (row['Hit'], row['Fam'])
        if fam_hit not in pairs and hit_fam not in pairs:
            unique_rows.append(row)
            pairs.add(fam_hit)

    hh_hits = pd.DataFrame(unique_rows)

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return hh_hits

def remove_tm(hh_hits, fams_to_export, tm_ids_file):
    with open(log_file, 'a') as file:
        file.write("Removing TM -")

    tm_ids = []
    with open(tm_ids_file, 'r') as file:
        for line in file:
            tm_ids.append(line.strip())

    with open(log_file, 'a') as file:
        file.write(f"{len(tm_ids)} fam(s)...")

    hh_hits        = hh_hits[~hh_hits['Fam'].isin(tm_ids) & ~hh_hits['Hit'].isin(tm_ids)]
    fams_to_export = [value for value in fams_to_export if value not in tm_ids]

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return hh_hits, fams_to_export

def read_fam_proteins_df(fam_proteins_file):
    with open(log_file, 'a') as file:
        file.write("Reading fam to proteins df...")

    column_names = ['Fam', 'Protein']
    fam_proteins = pd.read_csv(fam_proteins_file, sep='\t', header=None, names=column_names)
    fam_proteins[['Protein', 'Region']] = fam_proteins['Protein'].str.split('/', n=1, expand=True)

    with open(log_file, 'a') as file:
        file.write("Done\n")
    return fam_proteins

# Initiation Functions End ###

# Similarity Index Functions Start ###
def calculate_jaccard_index(set1, set2):
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    jaccard_index = len(intersection) / len(union)
    return jaccard_index

def get_fam_rep_subsets(fam1, fam2, fam_proteins):
    fam_rep = fam_to_rep_dict[fam1].split('/', 1)[0]
    
    fam_subset = fam_proteins[(fam_proteins['Fam'] == fam1) & (fam_proteins['Protein'] == fam_rep)][['Protein', 'Region']]
    hit_subset = fam_proteins[(fam_proteins['Fam'] == fam2) & (fam_proteins['Protein'] == fam_rep)][['Protein', 'Region']]

    return fam_subset, hit_subset

def create_aa_set(df, rep_length_dict):
    aa_set = set()

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        protein = row['Protein']
        region  = row['Region']
        
        # Split the region_str into individual regions
        regions = region.split('_')
        # Extract start and end numbers from the first and last region
        start_num = int(regions[0])
        end_num   = int(regions[-1])
        
        # Generate all combinations and add them to the set
        for num in range(start_num, end_num + 1):
            aa_set.add(f"{protein}:{num}")
            
    return aa_set

def calculate_aa_jaccard_index(fam, hit, fam_proteins, rep_length_dict):
    fam_subset, hit_subset = get_fam_rep_subsets(fam, hit, fam_proteins)
    id_to_remove = hit
    if hit_subset.empty: # if Fam rep protein name was not found in Hit
        fam_subset, hit_subset = get_fam_rep_subsets(hit, fam, fam_proteins)
        id_to_remove = fam
        if hit_subset.empty: # if Hit rep protein name was not found in Fam either, assume non redundant
            return -1, 0
    
    fam_subset['Region'] = fam_subset['Region'].fillna(f'1_{rep_length_dict[fam]}')
    hit_subset['Region'] = hit_subset['Region'].fillna(f'1_{rep_length_dict[fam]}')
    aa_set1              = create_aa_set(fam_subset, rep_length_dict)
    aa_set2              = create_aa_set(hit_subset, rep_length_dict)
    aa_jaccard_index     = calculate_jaccard_index(aa_set1, aa_set2)

    return aa_jaccard_index, id_to_remove

def remove_redundant(hh_hits, fam, redundant_fam_ids_file):
    hh_hits = hh_hits[(hh_hits['Fam'] != fam) & (hh_hits['Hit'] != fam)]

    with open(redundant_fam_ids_file, 'a') as file:
        file.write(f"{fam}\n")
    return hh_hits

def remove_pair(hh_hits, fam, hit):
    hh_hits = hh_hits[(hh_hits['Fam'] != fam) & (hh_hits['Hit'] != hit)]
    hh_hits = hh_hits[(hh_hits['Fam'] != hit) & (hh_hits['Hit'] != fam)]
    return hh_hits

def check_jaccard_similarity_remove_if_redundant(hh_hits, row, \
    fam_proteins, rep_length_dict, \
    fams_to_export, redundant_fam_ids_file, similarity_edgelist_file, \
    redundant_threshold, similarity_threshold):
    fam = row['Fam']
    hit = row['Hit']
    
    set_fam       = set(fam_proteins[fam_proteins['Fam'] == fam]['Protein'])
    set_hit       = set(fam_proteins[fam_proteins['Fam'] == hit]['Protein'])
    jaccard_index = calculate_jaccard_index(set_fam, set_hit)
    
    if (jaccard_index >= similarity_threshold): # 0.5
        with open(log_file, 'a') as file:
            file.write(f"{fam} vs {hit}: Initial Jaccard Index: {jaccard_index}...")
        aa_jaccard_index, id_to_remove = calculate_aa_jaccard_index(fam, hit, fam_proteins, rep_length_dict)

        if (aa_jaccard_index >= redundant_threshold): # 0.95
            with open(log_file, 'a') as file:
                file.write(f"AA Jaccard Index: {aa_jaccard_index}. Removing {id_to_remove}\n")
            fams_to_export.remove(id_to_remove)
            hh_hits = remove_redundant(hh_hits, id_to_remove, redundant_fam_ids_file)
        elif (aa_jaccard_index >= similarity_threshold): # 0.5
            with open(log_file, 'a') as file:
                file.write(f"AA Jaccard Index: {aa_jaccard_index}. Keeping similarity edge\n")
            with open(similarity_edgelist_file, 'a') as f: 
                f.write(f"{fam},{hit},{aa_jaccard_index}\n")
            hh_hits = remove_pair(hh_hits, fam, hit) # don't check same pair again
        else:
            with open(log_file, 'a') as file:
                file.write(f"AA Jaccard Index: {aa_jaccard_index}. Not similar enough\n")
    
    return hh_hits, fams_to_export

# Might still use this while comparing MGnifams to other resource families
# Don't call keep_unique(hh_hits) before this
def check_similarity_remove_if_redundant_hhblits_only(hh_hits, row, \
    fams_to_export, redundant_fam_ids_file, similarity_edgelist_file, \
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
                hh_hits = remove_redundant(hh_hits, hit, redundant_fam_ids_file)
            else:
                fams_to_export.remove(fam)
                hh_hits = remove_redundant(hh_hits, fam, redundant_fam_ids_file)
        else: # check Cols
            avg_cols1  = hh_hits[(hh_hits['Fam'] == fam) & (hh_hits['Hit'] == hit)]['Cols'].mean()
            avg_cols2  = hh_hits[(hh_hits['Fam'] == hit) & (hh_hits['Hit'] == fam)]['Cols'].mean()
            if (avg_cols1 >= avg_cols2):
                fams_to_export.remove(hit)
                hh_hits = remove_redundant(hh_hits, hit, redundant_fam_ids_file)
            else:
                fams_to_export.remove(fam)
                hh_hits = remove_redundant(hh_hits, fam, redundant_fam_ids_file)
    elif (avg_prob >= similarity_threshold):
        with open(similarity_edgelist_file, 'a') as f: 
            f.write(f"{fam},{hit},{avg_prob}\n")
        hh_hits = remove_pair(hh_hits, fam, hit)

    return hh_hits, fams_to_export

def write_non_redundant_fam_ids(fams_to_export, non_redundant_fam_ids_file):
    with open(non_redundant_fam_ids_file, 'w') as f:
        for id in fams_to_export:
            f.write(f"{id}\n")

def export_non_redundant_family_ids(hh_hits_file, fam_rep_mapping_file, \
    tm_ids_file, fam_proteins_file, rep_fa_file, \
    non_redundant_fam_ids_file, redundant_fam_ids_file, similarity_edgelist_file, log_f, \
    redundant_threshold=0.95, similarity_threshold=0.5):

    global rep_to_fam_dict, fam_to_rep_dict, log_file
    log_file                         = log_f
    rep_length_dict                  = create_fasta_to_length_dict(rep_fa_file)
    rep_to_fam_dict, fam_to_rep_dict = read_rep_to_fam_dicts(fam_rep_mapping_file)
    hh_hits                          = read_hh_hits(hh_hits_file)
    fams_to_export                   = hh_hits['Fam'].unique().tolist() # This must be done here, before removing self-hits (some fams might have only self-hits)
    hh_hits                          = map_and_remove_self(hh_hits)
    hh_hits                          = keep_unique_pairs(hh_hits)
    hh_hits, fams_to_export          = remove_tm(hh_hits, fams_to_export, tm_ids_file)
    fam_proteins                     = read_fam_proteins_df(fam_proteins_file)
    
    i = 0
    while len(hh_hits) > i:
        row = hh_hits.iloc[i]
        hh_hits, fams_to_export = check_jaccard_similarity_remove_if_redundant(hh_hits, row, \
            fam_proteins, rep_length_dict, \
            fams_to_export, redundant_fam_ids_file, similarity_edgelist_file, \
            redundant_threshold, similarity_threshold)
        i += 1

    write_non_redundant_fam_ids(fams_to_export, non_redundant_fam_ids_file)

if __name__ == "__main__":
    if len(sys.argv) != 10:
        print("Usage: python remove_redundant_and_tm.py <hh_hits> <fam_rep_mapping> <tm_ids_file> <fam_proteins_file> <rep_fa_file> \
            <non_redundant_fam_ids> <redundant_fam_ids> <similarity_edgelist> <log.txt>")
        sys.exit(1)

    initialize_outfiles(sys.argv[6:10]) # for the three output files: 6, 7, 8
    export_non_redundant_family_ids(sys.argv[1], sys.argv[2], \
        sys.argv[3], sys.argv[4], sys.argv[5], \
        sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
