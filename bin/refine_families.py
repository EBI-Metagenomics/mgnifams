import sys
import os
import subprocess
import shutil
import pickle
import pandas as pd
from Bio import SeqIO

import time # benchmarking, TODO remove

def parse_args():
    # Minimum 5 arguments (including the script name), maximum 7
    if not (5 <= len(sys.argv) <= 7):
        print("""
            Usage: python3 refine_families.py [Linclust file] [MGnifams FASTA file] [Minimum family members]
            [Output families file] [Optional 1: Clusters bookkeeping df .pkl file] [Optional 2: Updated mgnifams input dict .pkl file]
        """)
        sys.exit(1)

    globals().update({
        "linclust_input_file"    : sys.argv[1],
        "mgnifams_input_file"    : sys.argv[2],
        "minimum_family_members" : int(sys.argv[3]),
        "output_families_file"   : sys.argv[4]
    })

    # Optional arguments
    globals().update({"arg_clusters_bookkeeping_df_file": sys.argv[5] if len(sys.argv) >= 6 else None})
    globals().update({"arg_updated_mgnifams_input_file" : sys.argv[6] if len(sys.argv) >= 7 else None})        

def define_globals():
    globals().update({
        "log_file"                     : "log.txt",
        "clusters_bookkeeping_df_file" : "clusters_bookkeeping_df.pkl",
        "updated_mgnifams_input_file"  : "updated_mgnifams_input.pkl",
        "tmp_folder"                   : "tmp",
        "seed_msa_folder"              : "seed_msa",
        "align_msa_folder"             : "msa",
        "hmm_folder"                   : "hmm",
        "domtblout_folder"             : "dombtblout",
        "evalue_threshold"             : 0.001,
        "length_threshold"             : 0.8
    })

    for folder in [tmp_folder, seed_msa_folder, align_msa_folder, hmm_folder, domtblout_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    globals().update({
        "tmp_family_sequences_path"    : os.path.join(tmp_folder, 'family_sequences.fa'),
        "tmp_seed_msa_path"            : os.path.join(tmp_folder, 'seed_msa.fa'),
        "tmp_align_msa_path"           : os.path.join(tmp_folder, 'align_msa.fa'),
        "tmp_hmm_path"                 : os.path.join(tmp_folder, 'model.hmm'),
        "tmp_domtblout_path"           : os.path.join(tmp_folder, 'domtblout.txt'),
        "tmp_sequences_to_remove_path" : os.path.join(tmp_folder, 'sequences_to_remove.txt'),
        "seed_msa_path"                : os.path.join(seed_msa_folder, 'seed_msa.fa'),
        "align_msa_path"               : os.path.join(align_msa_folder, 'align_msa.fa'),
        "hmm_path"                     : os.path.join(hmm_folder, 'model.hmm'),
        "domtblout_path"               : os.path.join(domtblout_folder, 'domtblout.txt')
    })
    
def create_clusters_bookkeeping_df():
    start_time = time.time()

    df = pd.read_csv(linclust_input_file, sep='\t', header=None, names=['representative', 'member'])
    family_sizes = df.groupby('representative').size()
    clusters_bookkeeping_df = df.set_index('representative').join(family_sizes.rename('size'), on='representative')
    clusters_bookkeeping_df.to_pickle(clusters_bookkeeping_df_file)

    with open(log_file, 'w') as file:
        file.write("create_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def load_clusters_bookkeeping_df():
    start_time = time.time()

    clusters_bookkeeping_df = pd.read_pickle(arg_clusters_bookkeeping_df_file)

    with open(log_file, 'w') as file:
        file.write("load_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def create_mgnifams_input_dict():
    start_time = time.time()

    fasta_dict = {record.id: record for record in SeqIO.parse(mgnifams_input_file, "fasta")}
    with open(updated_mgnifams_input_file, 'wb') as file:
        pickle.dump(fasta_dict, file)

    with open(log_file, 'a') as file:
        file.write("create_mgnifams_input_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return fasta_dict

def load_mgnifams_input_dict():
    start_time = time.time()

    with open(arg_updated_mgnifams_input_file, 'rb') as file:
        fasta_dict = pickle.load(file)

    with open(log_file, 'a') as file:
        file.write("load_mgnifams_input_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return fasta_dict

def get_next_largest_family(clusters_bookkeeping_df):
    start_time = time.time()
    
    if clusters_bookkeeping_df.empty:
        return None, None

    largest_family_rep = clusters_bookkeeping_df['size'].idxmax()
    largest_family_members = clusters_bookkeeping_df.loc[largest_family_rep]['member']

    with open(log_file, 'a') as file:
        file.write("get_next_largest_family: ")
        file.write(str(time.time() - start_time) +"\n")
        file.write(f"S: {largest_family_rep}, s: {len(largest_family_members)}\n")

    return largest_family_rep, largest_family_members.tolist() if isinstance(largest_family_members, pd.Series) else [largest_family_members]

def write_fasta_file(members, fasta_dict, output_fasta, mode):
    start_time = time.time()
    sequences_to_write = [fasta_dict[member] for member in members if member in fasta_dict]

    with open(output_fasta, mode) as output_handle:
        SeqIO.write(sequences_to_write, output_handle, "fasta")

    with open(log_file, 'a') as file:
        file.write("write_fasta_file: ")
        file.write(str(time.time() - start_time) + "\n")

def run_msa(input_fasta, output_msa):
    start_time = time.time()
    
    mafft_command = ["mafft", "--quiet", "--retree", "2", "--maxiterate", "2", "--thread", "-1", input_fasta]
    with open(output_msa, "w") as output_handle:
        subprocess.run(mafft_command, stdout=output_handle)

    with open(log_file, 'a') as file:
        file.write("run_msa: ")
        file.write(str(time.time() - start_time) + "\n")

def run_hmmbuild(msa_file, output_hmm):
    start_time = time.time()                    

    hmmbuild_command = ["hmmbuild", "--informat", "afa", output_hmm, msa_file]
    subprocess.run(hmmbuild_command, stdout=subprocess.DEVNULL)
    
    with open(log_file, 'a') as file:
        file.write("run_hmmbuild: ")
        file.write(str(time.time() - start_time) + "\n")

def run_hmmsearch(hmm_file, fasta_file, output_recruitment):
    start_time = time.time()

    hmmsearch_command = ["hmmsearch", "--domtblout", output_recruitment, hmm_file, fasta_file]
    subprocess.run(hmmsearch_command, stdout=subprocess.DEVNULL)

    with open(log_file, 'a') as file:
        file.write("run_hmmsearch: ")
        file.write(str(time.time() - start_time) + "\n")

def filter_recruited(recruitment_file, evalue_threshold, length_threshold):
    start_time = time.time()

    filtered_sequences = []
    with open(recruitment_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split()
                evalue = float(columns[6])
                qlen = float(columns[5])
                hmm_from = int(columns[15])
                hmm_to = int(columns[16])
                hmm_length = hmm_to - hmm_from + 1
                if evalue < evalue_threshold and hmm_length >= length_threshold * qlen:
                    filtered_sequences.append(columns[0])
    
    with open(log_file, 'a') as file:
        file.write("filter_recruited: ")
        file.write(str(time.time() - start_time) + "\n")

    return filtered_sequences

def run_hmmalign(input_file, hmm_file, output_file):
    start_time = time.time()

    hmmalign_command = ["hmmalign", "--out", output_file, hmm_file, input_file]
    subprocess.run(hmmalign_command, stdout=subprocess.DEVNULL)

    with open(log_file, 'a') as file:
        file.write("run_hmmalign: ")
        file.write(str(time.time() - start_time) + "\n")

def run_esl_weight(input_file, output_file, threshold=0.80):
    start_time = time.time()

    esl_weight_command = ["esl-weight", "-o", output_file, "-t", str(threshold), input_file]
    subprocess.run(esl_weight_command, stdout=subprocess.DEVNULL)

    with open(log_file, 'a') as file:
        file.write("run_esl_weight: ")
        file.write(str(time.time() - start_time) + "\n")

def append_family_file(output_file, family_rep, family_members):
    start_time = time.time()

    lines = [f"{family_rep}\t{member}\n" for member in family_members]
    with open(output_file, 'a') as file:
        file.writelines(lines)

    with open(log_file, 'a') as file:
        file.write("append_family_file: ")
        file.write(str(time.time() - start_time) + "\n")
        file.write(f"E: {family_rep}, s: {len(family_members)}\n")

def save_iteration(clusters_bookkeeping_df, family_rep):
    start_time = time.time()

    os.rename(tmp_seed_msa_path, os.path.join(msa_folder, f"{family_rep}.msa"))
    os.rename(hmm_folder, os.path.join(hmm_folder, f"{family_rep}.hmm"))
    clusters_bookkeeping_df.to_pickle(tmp_bookkeeping_df_path)

    with open(log_file, 'a') as file:
        file.write("save_iteration: ")
        file.write(str(time.time() - start_time) + "\n")

def update_bookkeeping(filtered_seq_names, clusters_bookkeeping_df, fasta_dict):
    start_time = time.time()    

    with open(tmp_sequences_to_remove_path, "a") as to_remove_file:
        for seq_name in filtered_seq_names:
            to_remove_file.write(f"^>{seq_name}$\n")
            to_remove_file.write(f"^{fasta_dict[seq_name].seq}$\n")

            del fasta_dict[seq_name] # remove from dict

            representative = clusters_bookkeeping_df[clusters_bookkeeping_df['member'] == seq_name].index[0]
            clusters_bookkeeping_df.loc[clusters_bookkeeping_df.index.get_level_values('representative') == representative, 'size'] -= 1 # size -1 for family

    mask = ~clusters_bookkeeping_df['member'].isin(filtered_seq_names)
    clusters_bookkeeping_df = clusters_bookkeeping_df[mask].copy() # remove from df

    with open(log_file, 'a') as file:
        file.write("update_bookkeeping: ")
        file.write(str(time.time() - start_time) + "\n")

def remove_sequences_from_pool(to_remove_file, fasta_file):
    start_time = time.time()

    grep_command = ["grep", "-v", "-f", to_remove_file, fasta_file]
    with open(tmp_mgnifams_input_path, "w") as outfile:
        subprocess.run(grep_command, stdout=outfile)
    shutil.move(tmp_mgnifams_input_path, fasta_file)
    if os.path.exists(to_remove_file):
        os.remove(to_remove_file)

    with open(log_file, 'a') as file:
        file.write("remove_sequences_from_pool: ")
        file.write(str(time.time() - start_time) + "\n")

def main():
    parse_args()
    define_globals()

    clusters_bookkeeping_df = create_clusters_bookkeeping_df() if arg_clusters_bookkeeping_df_file is None else load_clusters_bookkeeping_df()
    fasta_dict = create_mgnifams_input_dict() if arg_updated_mgnifams_input_file is None else load_mgnifams_input_dict()

    exit() # TODO continue from here
    
    while True:
        family_rep, family_members = get_next_largest_family(clusters_bookkeeping_df)
        if not family_members or len(family_members) < minimum_family_members:
            break
        
        family_iteration = 0
        write_fasta_mode = "w"
        total_recruited_sequences = []
        while True:
            family_iteration += 1
            if (family_iteration > 3):
                break
            with open(log_file, 'a') as file:
                file.write(str(family_iteration) + "\n")

            write_fasta_file(family_members, fasta_dict, tmp_family_sequences_path, write_fasta_mode) # TODO check update
            run_msa(tmp_family_sequences_path, tmp_seed_msa_path)
            run_hmmbuild(tmp_seed_msa_path, tmp_hmm_path)
            run_hmmsearch(tmp_hmm_path, mgnifams_input_file, tmp_domtblout_path) # TODO update mgnifams_input_file used
            filtered_seq_names = filter_recruited(tmp_domtblout_path, evalue_threshold, length_threshold)
            recruited_sequences = set(filtered_seq_names) - set(family_members)
            if recruited_sequences:
                # TODO hmmalign
                run_hmmalign(input_file, hmm_file, output_file)
                break
                # TODO esl_weight, 0.8 until <2000

                total_recruited_sequences.extend(recruited_sequences)
                family_members.extend(recruited_sequences)
                write_fasta_mode = "a"
                continue
            else:
                break

        append_family_file(output_families_file, family_rep, family_members)
        save_iteration(clusters_bookkeeping_df, family_rep)
        update_bookkeeping(total_recruited_sequences, clusters_bookkeeping_df, fasta_dict)
        remove_sequences_from_pool(tmp_sequences_to_remove_path, mgnifams_input_file) # TODO update mgnifams_input_file used
        break # TODO export family stuff here

if __name__ == "__main__":
    main()
