import sys
import os
import subprocess
import shutil
import pickle
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO

import time # benchmarking

def parse_args():
    if not (len(sys.argv) == 4):
        print("Usage: python3 refine_families.py [Clusters bookkeeping df pkl file] [MGnifams FASTA file] [Number of minimum family members to check]")
        sys.exit(1)

    globals().update({
        "arg_clusters_bookkeeping_df_pkl_file" : sys.argv[1],
        "arg_mgnifams_dict_fasta_file"         : sys.argv[2],
        "arg_minimum_family_members"           : int(sys.argv[3])
    })   

def define_globals():
    globals().update({
        "output_families_file"                     : "refined_families.tsv",
        "log_file"                                 : "log.txt",
        "updated_clusters_bookkeeping_df_pkl_file" : "updated_clusters_bookkeeping_df.pkl",
        "updated_mgnifams_dict_fasta_file"         : "updated_mgnifams_dict.fa",
        "tmp_folder"                               : "tmp",
        "seed_msa_folder"                          : "seed_msa",
        "align_msa_folder"                         : "msa",
        "hmm_folder"                               : "hmm",
        "domtblout_folder"                         : "domtblout",
        "evalue_threshold"                         : 0.001,
        "length_threshold"                         : 0.8
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
        "tmp_intermediate_esl_path"    : os.path.join(tmp_folder, 'intermediate_esl.fa'),
        "tmp_esl_weight_path"          : os.path.join(tmp_folder, 'esl_weight.fa'),
        "tmp_sequences_to_remove_path" : os.path.join(tmp_folder, 'sequences_to_remove.txt')
    })

def load_clusters_bookkeeping_df():
    start_time = time.time()

    clusters_bookkeeping_df = pd.read_pickle(arg_clusters_bookkeeping_df_pkl_file)

    with open(log_file, 'w') as file:
        file.write("load_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def create_mgnifams_fasta_dict():
    start_time = time.time()

    mgnifams_fasta_dict = {record.id: record for record in SeqIO.parse(arg_mgnifams_dict_fasta_file, "fasta")}
    shutil.copy(arg_mgnifams_dict_fasta_file, updated_mgnifams_dict_fasta_file)

    with open(log_file, 'a') as file:
        file.write("create_mgnifams_fasta_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return mgnifams_fasta_dict

def get_next_largest_family(clusters_bookkeeping_df):
    start_time = time.time()
    
    if clusters_bookkeeping_df.empty:
        return None, None

    largest_family_rep = clusters_bookkeeping_df['size'].idxmax()
    largest_family_members = clusters_bookkeeping_df.loc[largest_family_rep]['member']
    if isinstance(largest_family_members, str):
        # If it's a single string, convert to a list with one element, else len counts the str characters
        largest_family_members = [largest_family_members]

    with open(log_file, 'a') as file:
        file.write("get_next_largest_family: ")
        file.write(str(time.time() - start_time) +"\n")
        file.write(f"S: {largest_family_rep}, s: {len(largest_family_members)}\n")

    return largest_family_rep, largest_family_members.tolist() if isinstance(largest_family_members, pd.Series) else [largest_family_members]

def write_family_fasta_file(members, mgnifams_fasta_dict, output_fasta):
    start_time = time.time()
    sequences_to_write = [mgnifams_fasta_dict[member] for member in members if member in mgnifams_fasta_dict]

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(sequences_to_write, output_handle, "fasta")

    with open(log_file, 'a') as file:
        file.write("write_family_fasta_file: ")
        file.write(str(time.time() - start_time) + "\n")

def run_msa(input_fasta, output_msa, family_size):
    start_time = time.time()
    
    if (family_size > 1):
        mafft_command = ["mafft", "--quiet", "--retree", "2", "--maxiterate", "2", "--thread", "-1", input_fasta]
    else:
        with open(log_file, 'a') as file:
            file.write("Running 1-sequence version of mafft. ")
        mafft_command = ["mafft", "--quiet", "--retree", "2", "--thread", "-1", input_fasta]

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

def mask_sequence_name(sequence_name, env_from, env_to, mgnifams_fasta_dict):
    masked_name = f"{sequence_name}/{env_from}_{env_to}"

    if sequence_name in mgnifams_fasta_dict:
        sub_sequence = str(mgnifams_fasta_dict[sequence_name].seq[env_from - 1:env_to])  # -1 because Python uses 0-based indexing
        masked_seq_record = SeqRecord(Seq(sub_sequence), id=masked_name, description="")

        with open(tmp_family_sequences_path, "a") as output_handle:
            SeqIO.write(masked_seq_record, output_handle, "fasta")

    return masked_name

def filter_recruited(recruitment_file, evalue_threshold, length_threshold, mgnifams_fasta_dict):
    start_time = time.time()

    os.remove(tmp_family_sequences_path)
    filtered_sequences = []
    with open(recruitment_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split()
                evalue = float(columns[6])
                qlen = float(columns[5])
                env_from = int(columns[19])
                env_to = int(columns[20])
                env_length = env_to - env_from + 1
                if evalue < evalue_threshold and env_length >= length_threshold * qlen:
                    sequence_name = columns[0]
                    tlen = float(columns[2])
                    if (env_length < tlen):
                        sequence_name = mask_sequence_name(sequence_name, env_from, env_to, mgnifams_fasta_dict)
                    else:
                        with open(tmp_family_sequences_path, "a") as output_handle:
                            SeqIO.write([mgnifams_fasta_dict[sequence_name]], output_handle, "fasta")

                    filtered_sequences.append(sequence_name)
    
    with open(log_file, 'a') as file:
        file.write("filter_recruited: ")
        file.write(str(time.time() - start_time) + "\n")

    return filtered_sequences

def run_hmmalign(input_file, hmm_file, output_file):
    start_time = time.time()

    hmmalign_command = ["hmmalign", "-o", output_file, "--amino", hmm_file, input_file]
    subprocess.run(hmmalign_command, stdout=subprocess.DEVNULL)

    with open(log_file, 'a') as file:
        file.write("run_hmmalign: ")
        file.write(str(time.time() - start_time) + "\n")

def get_number_of_remaining_sequences(output_file):
    with open(output_file, "r") as file:
        alignment = AlignIO.read(file, "stockholm")

    return len(alignment)

def run_esl_weight(input_file, output_file, threshold=0.8):
    start_time = time.time()

    shutil.copy(input_file, tmp_intermediate_esl_path)

    while True:
        esl_weight_command = ["esl-weight", "--amino", "-f", "--idf", str(threshold), "-o", output_file, tmp_intermediate_esl_path]
        subprocess.run(esl_weight_command, stdout=subprocess.DEVNULL)
        number_of_remaining_sequences = get_number_of_remaining_sequences(output_file)
        with open(log_file, 'a') as file:
            file.write("Remaining sequences: " + str(number_of_remaining_sequences) + "\n")
        if (number_of_remaining_sequences <= 2000):
            threshold -= 0.1
            break
        else:
            with open(log_file, 'a') as file:
                file.write("Rerunning run_esl_weight: ")
            shutil.copy(output_file, tmp_intermediate_esl_path)  

    with open(log_file, 'a') as file:
        file.write("run_esl_weight: ")
        file.write(str(time.time() - start_time) + "\n")

def filter_out_redundant(tmp_family_sequences_path, tmp_esl_weight_path):
    start_time = time.time()

    # Step 1: Read identifiers from tmp_esl_weight_path using AlignIO
    identifiers = set()
    with open(tmp_esl_weight_path, "r") as stockholm_file:
        alignment = AlignIO.read(stockholm_file, "stockholm")
        for record in alignment:
            identifiers.add(record.id)

    # Step 2: Read and filter sequences in tmp_family_sequences_path
    filtered_sequences = []
    kept_identifiers = [] 
    for record in SeqIO.parse(tmp_family_sequences_path, "fasta"):
        if record.id.split('/')[0] in identifiers:
            filtered_sequences.append(record)
            kept_identifiers.append(record.id)

    # Step 3: Write the filtered sequences back
    SeqIO.write(filtered_sequences, tmp_family_sequences_path, "fasta")

    with open(log_file, 'a') as file:
        file.write("filter_out_redundant: ")
        file.write(str(time.time() - start_time) + "\n")

    return kept_identifiers

def append_family_file(output_file, family_rep, family_members):
    start_time = time.time()

    lines = [f"{family_rep}\t{member}\n" for member in family_members]
    with open(output_file, 'a') as file:
        file.writelines(lines)

    with open(log_file, 'a') as file:
        file.write("append_family_file: ")
        file.write(str(time.time() - start_time) + "\n")
        file.write(f"E: {family_rep}, s: {len(family_members)}\n")

def move_produced_models(family_rep, size):
    shutil.move(tmp_seed_msa_path, os.path.join(seed_msa_folder, f'{family_rep}_{size}.fa'))
    shutil.move(tmp_align_msa_path, os.path.join(align_msa_folder, f'{family_rep}_{size}.fa'))
    shutil.move(tmp_hmm_path, os.path.join(hmm_folder, f'{family_rep}_{size}.hmm'))
    shutil.move(tmp_domtblout_path, os.path.join(domtblout_folder, f'{family_rep}_{size}.domtblout'))

def update_clusters_bookkeeping_df(clusters_bookkeeping_df, family_members):
    start_time = time.time()    

    sequences_to_remove = []
    split_members = [member.split('/')[0] for member in family_members]
    # Match these sequences to their family representatives and keep unique
    unique_reps = clusters_bookkeeping_df[clusters_bookkeeping_df['member'].isin(split_members)].index.unique()
    # Get all members of these unique reps into sequences_to_remove
    for rep in unique_reps:
        members = clusters_bookkeeping_df.loc[rep, 'member']
        if isinstance(members, pd.Series):
            sequences_to_remove.extend(members.tolist())
        else:
            sequences_to_remove.append(members)  
    # Remove all lines having these family_reps from the clusters_bookkeeping_df
    clusters_bookkeeping_df = clusters_bookkeeping_df[~clusters_bookkeeping_df.index.isin(unique_reps)]
    
    clusters_bookkeeping_df.to_pickle(updated_clusters_bookkeeping_df_pkl_file)

    with open(log_file, 'a') as file:
        file.write("update_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df, sequences_to_remove

def update_mgnifams_fasta_dict(mgnifams_fasta_dict, sequences_to_remove):
    start_time = time.time()
    
    # Filter out the sequences to remove using dictionary comprehension
    mgnifams_fasta_dict = {key: value for key, value in mgnifams_fasta_dict.items() if key not in sequences_to_remove}

    with open(updated_mgnifams_dict_fasta_file, 'w') as out_file:
        SeqIO.write(mgnifams_fasta_dict.values(), out_file, 'fasta')

    with open(log_file, 'a') as file:
        file.write("update_mgnifams_fasta_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return mgnifams_fasta_dict

def remove_tmp_files(folder_path):
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isfile(item_path):
            os.remove(item_path)

def main():
    parse_args()
    define_globals()

    clusters_bookkeeping_df = load_clusters_bookkeeping_df()
    mgnifams_fasta_dict = create_mgnifams_fasta_dict()
    while True:
        family_rep, family_members = get_next_largest_family(clusters_bookkeeping_df)
        if not family_members or len(family_members) < arg_minimum_family_members:
            with open(log_file, 'a') as file:
                file.write("Exiting all...")
            break
        
        write_family_fasta_file(family_members, mgnifams_fasta_dict, tmp_family_sequences_path)
        total_checked_sequences = family_members
        exit_flag = False
        family_iteration = 0
        while True:
            family_iteration += 1
            if (family_iteration > 3):
                exit_flag = True

            with open(log_file, 'a') as file:
                if (exit_flag):
                    file.write("Exiting family...: ")
                file.write(str(family_iteration) + "\n")
            
            run_msa(tmp_family_sequences_path, tmp_seed_msa_path, len(family_members))
            run_hmmbuild(tmp_seed_msa_path, tmp_hmm_path)
            run_hmmsearch(tmp_hmm_path, updated_mgnifams_dict_fasta_file, tmp_domtblout_path)
            filtered_seq_names = filter_recruited(tmp_domtblout_path, evalue_threshold, length_threshold, mgnifams_fasta_dict) # also writes in tmp_family_sequences_path
            run_hmmalign(tmp_family_sequences_path, tmp_hmm_path, tmp_align_msa_path)
            
            recruited_sequences = set(filtered_seq_names) - set(total_checked_sequences)
            if not recruited_sequences:
                exit_flag = True
            if (exit_flag):
                break

            total_checked_sequences.extend(filtered_seq_names)
            total_checked_sequences = list(set(total_checked_sequences))
            run_esl_weight(tmp_align_msa_path, tmp_esl_weight_path)
            family_members = filter_out_redundant(tmp_family_sequences_path, tmp_esl_weight_path) # also writes in tmp_family_sequences_path

        # Exiting family loop
        append_family_file(output_families_file, family_rep, family_members)
        move_produced_models(family_rep, len(family_members))
        clusters_bookkeeping_df, sequences_to_remove = update_clusters_bookkeeping_df(clusters_bookkeeping_df, family_members)
        mgnifams_fasta_dict = update_mgnifams_fasta_dict(mgnifams_fasta_dict, sequences_to_remove)
        remove_tmp_files(tmp_folder)

    # End of all families

if __name__ == "__main__":
    main()
