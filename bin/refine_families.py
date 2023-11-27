import sys
import os
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO

import time # benchmarking, TODO remove

def parse_args():
        if len(sys.argv) != 5:
            print("Usage: python3 refine_families.py [Families TSV] [FASTA file] [Minimum number of family members] [Output file]")
            sys.exit(1)

        global families_tsv 
        families_tsv = sys.argv[1]
        global fasta_file
        fasta_file = sys.argv[2]
        global minimum_members
        minimum_members = int(sys.argv[3])
        global output_file
        output_file = sys.argv[4]

def define_globals():
    global log_file
    log_file = "log.txt"
    global tmp_folder
    tmp_folder = "tmp"
    global msa_folder
    msa_folder = "msa"
    global hmm_folder
    hmm_folder = "hmm"
    for folder in [tmp_folder, msa_folder, hmm_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
        if not os.path.exists(tmp_folder):
            os.makedirs(tmp_folder)
    global family_sequences_path
    family_sequences_path = os.path.join(tmp_folder, 'family_sequences.fasta')
    global family_alignment_path
    family_alignment_path = os.path.join(tmp_folder, 'family_alignment.msa')
    global family_model_path
    family_model_path = os.path.join(tmp_folder, 'family_model.hmm')
    global recruited_sequences_path
    recruited_sequences_path = os.path.join(tmp_folder, 'domtblout.txt')

    global evalue_threshold
    evalue_threshold = 0.001
    global length_threshold
    length_threshold = 0.8
    global recruited_sequences
    recruited_sequences = []
    
def create_family_dataframe(families_tsv):
    start_time = time.time()

    df = pd.read_csv(families_tsv, sep='\t', header=None, names=['representative', 'member'])
    family_sizes = df.groupby('representative').size()
    bookkeeping_df = df.set_index('representative').join(family_sizes.rename('size'), on='representative')

    with open(log_file, 'w') as file:
        file.write("create_family_dataframe: ")
        file.write(str(time.time() - start_time) +"\n")

    return bookkeeping_df

# TODO load state alternative
def load_mgnifams_input():
    start_time = time.time()
    with open(log_file, 'a') as file:
        file.write("load_mgnifams_input: ")

    fasta_dict = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

    with open(log_file, 'a') as file:
        file.write(str(time.time() - start_time) +"\n")

    return fasta_dict

def get_next_largest_family(bookkeeping_df):
    start_time = time.time()
    
    if bookkeeping_df.empty:
        return None, None

    largest_family_rep = bookkeeping_df['size'].idxmax()
    largest_family_members = bookkeeping_df.loc[largest_family_rep]['member']

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

def append_family_file(output_file, family_rep, family_members):
    start_time = time.time()

    lines = [f"{family_rep}\t{member}\n" for member in family_members]
    with open(output_file, 'a') as file:
        file.writelines(lines)

    with open(log_file, 'a') as file:
        file.write("append_family_file: ")
        file.write(str(time.time() - start_time) + "\n")
        file.write(f"E: {family_rep}, s: {len(family_members)}\n")

def save_iteration(bookkeeping_df):
    start_time = time.time()

    os.rename(family_alignment_path, os.path.join(msa_folder, f"{family_rep}.msa"))
    os.rename(family_model_path, os.path.join(hmm_folder, f"{family_rep}.hmm"))
    bookkeeping_df.to_pickle('tmp/bookkeeping_df.pkl')

    with open(log_file, 'a') as file:
        file.write("save_iteration: ")
        file.write(str(time.time() - start_time) + "\n")

def update_bookkeeping(filtered_seq_names, bookkeeping_df, fasta_dict):
    start_time = time.time()    

    with open("tmp/sequences_to_remove.txt", "a") as to_remove_file:
        for seq_name in filtered_seq_names:
            to_remove_file.write(f"^>{seq_name}$\n")
            to_remove_file.write(f"^{fasta_dict[seq_name].seq}$\n")

            del fasta_dict[seq_name] # remove from dict

            representative = bookkeeping_df[bookkeeping_df['member'] == seq_name].index[0]
            bookkeeping_df.loc[bookkeeping_df.index.get_level_values('representative') == representative, 'size'] -= 1 # size -1 for family

    mask = ~bookkeeping_df['member'].isin(filtered_seq_names)
    bookkeeping_df = bookkeeping_df[mask].copy() # remove from df

    with open(log_file, 'a') as file:
        file.write("update_bookkeeping: ")
        file.write(str(time.time() - start_time) + "\n")

def remove_sequences_from_pool(to_remove_file, fasta_file):
    start_time = time.time()

    grep_command = ["grep", "-v", "-f", to_remove_file, fasta_file]
    with open("tmp/mgnifams_input.fa", "w") as outfile:
        subprocess.run(grep_command, stdout=outfile)
    shutil.move("tmp/mgnifams_input.fa", fasta_file)
    if os.path.exists(to_remove_file):
        os.remove(to_remove_file)

    with open(log_file, 'a') as file:
        file.write("remove_sequences_from_pool: ")
        file.write(str(time.time() - start_time) + "\n")

def main():
    parse_args()
    define_globals()

    bookkeeping_df = create_family_dataframe(families_tsv)
    fasta_dict = load_mgnifams_input()
    
    while True:
        family_rep, family_members = get_next_largest_family(bookkeeping_df)
        if not family_members or len(family_members) < minimum_members:
            break
        
        family_iteration = 0
        write_fasta_mode = "w"
        total_recruited_sequences = []
        while True:
            family_iteration += 1
            with open(log_file, 'a') as file:
                file.write(str(family_iteration) + "\n")

            write_fasta_file(family_members, fasta_dict, family_sequences_path, write_fasta_mode)
            run_msa(family_sequences_path, family_alignment_path)
            # TODO esl-weight here if >1000 sequences in MSA?
            run_hmmbuild(family_alignment_path, family_model_path)
            run_hmmsearch(family_model_path, fasta_file, recruited_sequences_path)
            filtered_seq_names = filter_recruited(recruited_sequences_path, evalue_threshold, length_threshold)

            recruited_sequences = set(filtered_seq_names) - set(family_members)
            if recruited_sequences:
                total_recruited_sequences.extend(recruited_sequences)
                family_members.extend(recruited_sequences)
                write_fasta_mode = "a"
                continue
            else:
                break # TODO continue from here
                append_family_file(output_file, family_rep, family_members)
                save_iteration(bookkeeping_df)
                update_bookkeeping(total_recruited_sequences, bookkeeping_df, fasta_dict)
                remove_sequences_from_pool("tmp/sequences_to_remove.txt", fasta_file)
                break

        break

if __name__ == "__main__":
    main()
