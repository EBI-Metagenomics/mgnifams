import sys
import pandas as pd
import csv

def read_hh_hits(hh_hits_file):
    column_names = ["Fam", "Hit", "Prob", "E-value", "P-value", "Score", "SS", "Cols", "Query HMM", "Template HMM", ""]
    hh_hits      = pd.read_csv(hh_hits_file, delimiter='\t', header=None, names=column_names, skiprows=1)

    return hh_hits

def read_rep_to_fam_dicts(fam_rep_mapping_file):
    rep_to_fam_dict = {}
    fam_to_rep_dict = {}
    
    with open(fam_rep_mapping_file, mode='r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            rep = row[1]
            fam = row[0]
            rep_to_fam_dict[rep] = fam
            fam_to_rep_dict[fam] = rep

    return rep_to_fam_dict, fam_to_rep_dict

def map_and_remove_self(hh_hits):
    hh_hits['Hit'] = hh_hits['Hit'].map(rep_to_fam_dict) # Do the mapping on the df
    hh_hits = hh_hits[hh_hits['Fam'] != hh_hits['Hit']]
    return hh_hits

def remove_tm(hh_hits, fams_to_export, tm_ids_file):
    tm_ids = []
    with open(tm_ids_file, 'r') as file:
        for line in file:
            tm_ids.append(line.strip())

    hh_hits        = hh_hits[~hh_hits['Fam'].isin(tm_ids) & ~hh_hits['Hit'].isin(tm_ids)]
    fams_to_export = [value for value in fams_to_export if value not in tm_ids]
    return hh_hits, fams_to_export

def read_fam_proteins(fam_proteins_file):
    column_names = ['Fam', 'Protein']
    fam_proteins = pd.read_csv(fam_proteins_file, sep='\t', header=None, names=column_names)
    fam_proteins[['Protein', 'Region']] = fam_proteins['Protein'].str.split('/', n=1, expand=True)
    return fam_proteins

def remove_redundant(hh_hits, fam):
    hh_hits = hh_hits[(hh_hits['Fam'] != fam) & (hh_hits['Hit'] != fam)]
    return hh_hits

def remove_pair(hh_hits, fam, hit):
    hh_hits = hh_hits[(hh_hits['Fam'] != fam) & (hh_hits['Hit'] != hit)]
    hh_hits = hh_hits[(hh_hits['Fam'] != hit) & (hh_hits['Hit'] != fam)]
    return hh_hits

def find_max_lengths(df):
    max_lengths = {}

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        fam = row['Fam']
        hit = row['Hit']

        # Extract and process Query HMM for Fam
        if fam not in max_lengths:
            max_lengths[fam] = 0
        fam_query_hmm = row['Query HMM']
        if pd.notna(fam_query_hmm):
            fam_length = int(fam_query_hmm.split('-')[1])
            max_lengths[fam] = max(max_lengths[fam], fam_length)

        # Extract and process Template HMM for Hit
        if hit not in max_lengths:
            max_lengths[hit] = 0
        hit_template_hmm = row['Template HMM']
        if pd.notna(hit_template_hmm):
            hit_length = int(hit_template_hmm.split('-')[1])
            max_lengths[hit] = max(max_lengths[hit], hit_length)

    return max_lengths

def keep_unique_pairs(hh_hits):
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

    return hh_hits

def calculate_jaccard_index(set_fam, set_hit):
    intersection = set_fam.intersection(set_hit)
    union = set_fam.union(set_hit)
    jaccard_index = len(intersection) / len(union)
    return jaccard_index

def filter_protein_region_df(fam1, fam2, fam_proteins):
    fam_rep_whole = fam_to_rep_dict[fam1]
    if '/' in fam_rep_whole:
        fam_rep, region = fam_rep_whole.split('/', 1)
    else:
        fam_rep = fam_rep_whole
        region = None

    filtered_df = fam_proteins[(fam_proteins['Fam'] == fam2) & (fam_proteins['Protein'] == fam_rep)][['Protein', 'Region']]
    return filtered_df

def calculate_whole_region(protein, max_hmm_lengths):
    if ('_' in protein): # whole slice
        parts      = protein.split('_')
        region_str = f'{parts[1]}_{parts[2]}'
    else: # whole protein
        fam = rep_to_fam_dict[protein]
        region_str = f'1_{max_hmm_lengths[fam]}'

    return region_str

def create_aa_set(df, max_hmm_lengths):
    aa_set = set()

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        protein    = row['Protein']
        region_str = row['Region']
        
        if (region_str is None):
            region_str = calculate_whole_region(protein, max_hmm_lengths)
        
        # Split the region_str into individual regions
        regions = region_str.split('_')
        
        # Extract start and end numbers from the first and last region
        start_num = int(regions[0])
        end_num = int(regions[-1])
        
        # Generate all combinations and add them to the set
        for num in range(start_num, end_num + 1):
            aa_set.add(f"{protein}:{num}")
            
    return aa_set

def calculate_name_and_region(fam_rep_whole, max_hmm_lengths):
    fam_rep = fam_rep_whole
    if '/' in fam_rep_whole: # part of whole or slice
        fam_rep, region = fam_rep_whole.split('/', 1)
    elif ('_' in fam_rep_whole): # whole slice
        parts  = fam_rep_whole.split('_')
        region = f'{parts[1]}_{parts[2]}'
    else: # whole protein
        fam = rep_to_fam_dict[fam_rep_whole]
        region = f'1_{max_hmm_lengths[fam]}'
    return fam_rep, region

def calculate_aa_jaccard_index(fam, hit, fam_proteins, max_hmm_lengths):
    print("### calculate_aa_jaccard_index ###")
    print(f"{fam} vs {hit}")
    filtered_hit_df = filter_protein_region_df(fam, hit, fam_proteins)
    fam_rep_whole   = fam_to_rep_dict[fam]
    id_to_remove    = hit
    if filtered_hit_df.empty: # if Fam rep protein name was not found in Hit
        filtered_hit_df = filter_protein_region_df(hit, fam, fam_proteins)
        fam_rep_whole   = fam_to_rep_dict[hit]
        id_to_remove    = fam
    
        if filtered_hit_df.empty: # if Hit rep protein name was not found in Fam either, assume non redundant
            return 0, 0
    
    aa_set1          = create_aa_set(filtered_hit_df, max_hmm_lengths)
    fam_rep, region  = calculate_name_and_region(fam_rep_whole, max_hmm_lengths)
    df               = pd.DataFrame({'Protein': [fam_rep], 'Region': [region]})
    aa_set2          = create_aa_set(df, max_hmm_lengths)
    aa_jaccard_index = calculate_jaccard_index(aa_set1, aa_set2)
    print(f'### AA similarity: {aa_jaccard_index}, id to remove: {id_to_remove} ###')
    return aa_jaccard_index, id_to_remove

def check_jaccard_similarity_remove_if_redundant(hh_hits, row, \
    fam_proteins, max_hmm_lengths, \
    fams_to_export, similarity_edgelist_file, \
    redundant_threshold, similarity_threshold):
    fam = row['Fam']
    hit = row['Hit']
    
    set_fam       = set(fam_proteins[fam_proteins['Fam'] == fam]['Protein'])
    set_hit       = set(fam_proteins[fam_proteins['Fam'] == hit]['Protein'])
    jaccard_index = calculate_jaccard_index(set_fam, set_hit)
    
    if (jaccard_index >= similarity_threshold): # 0.5
        aa_jaccard_index, id_to_remove = calculate_aa_jaccard_index(fam, hit, fam_proteins, max_hmm_lengths)
        if (aa_jaccard_index >= redundant_threshold): # 0.95
            print(f"{fam} {hit} -> Removing {id_to_remove}")
            fams_to_export.remove(id_to_remove)
            hh_hits = remove_redundant(hh_hits, id_to_remove)
        elif (aa_jaccard_index >= similarity_threshold): # 0.5
            with open(similarity_edgelist_file, 'a') as f: 
                f.write(f"{fam},{hit},{aa_jaccard_index}\n")
            hh_hits = remove_pair(hh_hits, fam, hit)
    
    return hh_hits, fams_to_export

# Might still use this while comparing MGnifams to other resource families
# Don't call keep_unique(hh_hits) before this
def check_similarity_remove_if_redundant_hhblits_only(hh_hits, row, \
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
    tm_ids_file, fam_proteins_file, \
    non_redundant_fam_ids_file, similarity_edgelist_file, \
    redundant_threshold=0.95, similarity_threshold=0.5):

    global rep_to_fam_dict, fam_to_rep_dict

    hh_hits                          = read_hh_hits(hh_hits_file)
    fams_to_export                   = hh_hits['Fam'].unique().tolist() # This must be done here, before removing self-hits (some fams might have only self-hits)
    rep_to_fam_dict, fam_to_rep_dict = read_rep_to_fam_dicts(fam_rep_mapping_file)
    hh_hits                          = map_and_remove_self(hh_hits)
    hh_hits, fams_to_export          = remove_tm(hh_hits, fams_to_export, tm_ids_file)
    fam_proteins                     = read_fam_proteins(fam_proteins_file)
    max_hmm_lengths                  = find_max_lengths(hh_hits)
    hh_hits                          = keep_unique_pairs(hh_hits)
    
    i = 0
    while len(hh_hits) > i:
        row = hh_hits.iloc[i]
        hh_hits, fams_to_export = check_jaccard_similarity_remove_if_redundant(hh_hits, row, \
            fam_proteins, max_hmm_lengths, \
            fams_to_export, similarity_edgelist_file, \
            redundant_threshold, similarity_threshold)
        i += 1

    write_non_redundant_fam_ids(fams_to_export)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python remove_redundant_and_tm.py <hh_hits> <fam_rep_mapping> <tm_ids_file> <fam_proteins_file>")
        sys.exit(1)

    non_redundant_fam_ids_file = "non_redundant_fam_ids.txt"
    similarity_edgelist_file   = "similarity_edgelist.csv"
    with open("similarity_edgelist.csv", "w") as file:
        pass # init empty
    export_non_redundant_family_ids(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], non_redundant_fam_ids_file, similarity_edgelist_file)
