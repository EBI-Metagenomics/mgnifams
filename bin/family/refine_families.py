import sys
import os
import subprocess
import shutil
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
import time

def parse_args():
    global arg_clusters_bookkeeping_df_pkl_file, arg_refined_families_tsv_file, \
        arg_mgnifams_input_fasta_file, arg_discarded_clusters_file, \
        arg_converged_families_file, arg_family_metadata_file, \
        arg_starting_num_sequences, arg_minimum_family_members, \
        arg_iteration
        
    if not (len(sys.argv) == 10):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_clusters_bookkeeping_df_pkl_file = sys.argv[1]
    arg_refined_families_tsv_file        = sys.argv[2]
    arg_mgnifams_input_fasta_file        = sys.argv[3]
    arg_discarded_clusters_file          = sys.argv[4]
    arg_converged_families_file          = sys.argv[5]
    arg_family_metadata_file             = sys.argv[6]
    arg_starting_num_sequences           = sys.argv[7]
    arg_minimum_family_members           = int(sys.argv[8])
    arg_iteration                        = int(sys.argv[9])

def define_globals():
    global log_file, updated_refined_families_tsv_file, \
        updated_mgnifams_input_fasta_file, updated_discarded_clusters_file, \
        updated_converged_families_file, updated_family_metadata_file, \
        tmp_folder, seed_msa_folder, \
        align_msa_folder, hmm_folder, \
        domtblout_folder, rf_folder, \
        evalue_threshold, length_threshold, \
        tmp_family_sequences_path, tmp_seed_msa_path, \
        tmp_seed_msa_sto_path, tmp_align_msa_path, \
        tmp_hmm_path, tmp_domtblout_path, \
        tmp_intermediate_esl_path, tmp_esl_weight_path, \
        tmp_sequences_to_remove_path, tmp_rf_path

    log_file                          = "log.txt"
    updated_refined_families_tsv_file = "updated_refined_families.tsv"
    updated_mgnifams_input_fasta_file = "updated_mgnifams_input.fa"
    updated_discarded_clusters_file   = "updated_discarded_clusters.txt"
    updated_converged_families_file   = "updated_converged_families.txt"
    updated_family_metadata_file      = "updated_family_metadata.csv"
    tmp_folder                        = "tmp"
    seed_msa_folder                   = "seed_msa_sto"
    align_msa_folder                  = "msa_sto"
    hmm_folder                        = "hmm"
    domtblout_folder                  = "domtblout"
    rf_folder                         = "rf"
    evalue_threshold                  = 0.001
    length_threshold                  = 0.8

    for folder in [tmp_folder, seed_msa_folder, align_msa_folder, hmm_folder, domtblout_folder, rf_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    tmp_family_sequences_path    = os.path.join(tmp_folder, 'family_sequences.fa')
    tmp_seed_msa_path            = os.path.join(tmp_folder, 'seed_msa.fa')
    tmp_seed_msa_sto_path        = os.path.join(tmp_folder, 'seed_msa_sto.sto')
    tmp_align_msa_path           = os.path.join(tmp_folder, 'align_msa.sto')
    tmp_hmm_path                 = os.path.join(tmp_folder, 'model.hmm')
    tmp_domtblout_path           = os.path.join(tmp_folder, 'domtblout.txt')
    tmp_intermediate_esl_path    = os.path.join(tmp_folder, 'intermediate_esl.fa')
    tmp_esl_weight_path          = os.path.join(tmp_folder, 'esl_weight.fa')
    tmp_sequences_to_remove_path = os.path.join(tmp_folder, 'sequences_to_remove.txt')
    tmp_rf_path                  = os.path.join(tmp_folder, 'rf.txt')

def copy_updated_input_files():

    start_time = time.time()

    shutil.copy(arg_refined_families_tsv_file, updated_refined_families_tsv_file)
    shutil.copy(arg_mgnifams_input_fasta_file, updated_mgnifams_input_fasta_file)
    shutil.copy(arg_discarded_clusters_file, updated_discarded_clusters_file)
    shutil.copy(arg_converged_families_file, updated_converged_families_file)
    shutil.copy(arg_family_metadata_file, updated_family_metadata_file)

    with open(log_file, 'w') as file:
        file.write("copy_updated_input_files: ")
        file.write(str(time.time() - start_time) + "\n")

def load_clusters_bookkeeping_df():
    start_time = time.time()

    clusters_bookkeeping_df = pd.read_pickle(arg_clusters_bookkeeping_df_pkl_file)

    with open(log_file, 'a') as file:
        file.write("load_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def exclude_existing_families(clusters_bookkeeping_df):
    start_time = time.time()

    if os.path.getsize(arg_refined_families_tsv_file) > 0:
        refined_families_df = pd.read_csv(arg_refined_families_tsv_file, delimiter='\t', header=None, dtype=str)
        members_to_remove = refined_families_df.iloc[:, 1].apply(lambda x: x.split('/')[0]).unique().tolist()
        clusters_to_remove = clusters_bookkeeping_df[clusters_bookkeeping_df['member'].isin(members_to_remove)].index.unique().tolist()
        # Modify the DataFrame in place
        clusters_bookkeeping_df.drop(index=clusters_to_remove, inplace=True)

    with open(log_file, 'a') as file:
        file.write("exclude_existing_families: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def exclude_discarded_clusters(clusters_bookkeeping_df):
    start_time = time.time()

    if os.path.getsize(arg_discarded_clusters_file) > 0:
        discarded_clusters_df = pd.read_csv(arg_discarded_clusters_file, header=None, dtype=str)
        clusters_to_remove = clusters_bookkeeping_df[clusters_bookkeeping_df['member'].isin(discarded_clusters_df[0])].index.unique().tolist()
        # Modify the DataFrame in place
        clusters_bookkeeping_df.drop(index=clusters_to_remove, inplace=True)

    with open(log_file, 'a') as file:
        file.write("exclude_discarded_clusters: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def create_mgnifams_fasta_dict():
    start_time = time.time()

    mgnifams_fasta_dict = {record.id: record for record in SeqIO.parse(arg_mgnifams_input_fasta_file, "fasta")}

    with open(log_file, 'a') as file:
        file.write("create_mgnifams_fasta_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return mgnifams_fasta_dict

def get_next_largest_family(clusters_bookkeeping_df):
    start_time = time.time()
    
    if clusters_bookkeeping_df.empty:
        return None, None

    largest_family_rep = clusters_bookkeeping_df['size'].idxmax()
    largest_family_members = clusters_bookkeeping_df.loc[largest_family_rep, 'member']
    # Ensure that largest_family_members is a list
    if not isinstance(largest_family_members, list):
        largest_family_members = [largest_family_members] if isinstance(largest_family_members, str) else largest_family_members.tolist()

    with open(log_file, 'a') as file:
        file.write("\nget_next_largest_family: ")
        file.write(str(time.time() - start_time) + "\n")
        file.write(f"S: {largest_family_rep}, s: {len(largest_family_members)}\n")

    return largest_family_rep, largest_family_members

def write_fasta_sequences(sequences, file_path, mode):
    with open(file_path, mode) as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def write_family_fasta_file(members, mgnifams_fasta_dict):
    sequences_to_write = [mgnifams_fasta_dict[member] for member in members if member in mgnifams_fasta_dict]
    write_fasta_sequences(sequences_to_write, tmp_family_sequences_path, "w")

def run_msa(family_size):
    start_time = time.time()
    
    if (family_size > 1):
        mafft_command = ["mafft", "--quiet", "--retree", "2", "--maxiterate", "2", "--thread", "-1", tmp_family_sequences_path]
    else:
        with open(log_file, 'a') as file:
            file.write("Running 1-sequence version of mafft. ")
        mafft_command = ["mafft", "--quiet", "--retree", "2", "--thread", "-1", tmp_family_sequences_path]

    with open(tmp_seed_msa_path, "w") as output_handle:
        subprocess.run(mafft_command, stdout=output_handle)

    with open(log_file, 'a') as file:
        file.write("run_msa: ")
        file.write(str(time.time() - start_time) + "\n")

def read_fasta_to_matrix(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    max_length = max(len(record.seq) for record in records)
    matrix = np.zeros((len(records), max_length), dtype=np.dtype('U1'))
    original_names = []

    for i, record in enumerate(records):
        original_names.append(record.id)
        matrix[i, :len(record.seq)] = list(str(record.seq))

    return matrix, original_names

def calculate_trim_positions(sequence_matrix, occupancy_threshold):
    numeric_matrix = np.where(sequence_matrix == '-', 0, 1)
    num_rows = numeric_matrix.shape[0]
    column_sums = np.sum(numeric_matrix, axis=0)
    column_sums_percentage = column_sums / num_rows
    start_position = np.argmax(column_sums_percentage > occupancy_threshold)
    end_position = len(column_sums_percentage) - np.argmax(column_sums_percentage[::-1] > occupancy_threshold) - 1

    return start_position, end_position

def write_trimmed_sequences(sequence_matrix_trimmed, original_sequence_names):
    trimmed_records = []
    for i, sequence in enumerate(sequence_matrix_trimmed):
        trimmed_sequence = ''.join(map(str, sequence))
        original_name = original_sequence_names[i]
        trimmed_record = SeqIO.SeqRecord(Seq(trimmed_sequence), id=original_name, description="")
        trimmed_records.append(trimmed_record)

    with open(tmp_seed_msa_path, "w") as output_fasta:
        SeqIO.write(trimmed_records, output_fasta, "fasta")
        
def trim_seed_msa(occupancy_threshold=0.5):
    start_time = time.time()

    sequence_matrix, original_sequence_names = read_fasta_to_matrix(tmp_seed_msa_path)
    start_position, end_position = calculate_trim_positions(sequence_matrix, occupancy_threshold)
    sequence_matrix_trimmed = sequence_matrix[:, start_position:end_position+1]
    write_trimmed_sequences(sequence_matrix_trimmed, original_sequence_names)
    
    with open(log_file, 'a') as file:
        file.write("trim_seed_msa: ")
        file.write(str(time.time() - start_time) + "\n")

    return original_sequence_names

def run_hmmbuild(msa_file, extra_args):
    start_time = time.time()                    

    hmmbuild_command = ["hmmbuild"] + extra_args + [tmp_hmm_path, msa_file]
    subprocess.run(hmmbuild_command, stdout=subprocess.DEVNULL)
    
    with open(log_file, 'a') as file:
        file.write("run_hmmbuild: ")
        file.write(str(time.time() - start_time) + "\n")

def run_hmmsearch():
    start_time = time.time()
    
    hmmsearch_command = ["hmmsearch", "-Z", arg_starting_num_sequences, "--domtblout", tmp_domtblout_path, tmp_hmm_path, updated_mgnifams_input_fasta_file]
    subprocess.run(hmmsearch_command, stdout=subprocess.DEVNULL)

    with open(log_file, 'a') as file:
        file.write("run_hmmsearch: ")
        file.write(str(time.time() - start_time) + "\n")

def extract_sequence_names_from_domtblout():
    sequence_names = []

    with open(tmp_domtblout_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue
            columns = line.split()
            sequence_name = columns[0]
            sequence_names.append(sequence_name)

    return sequence_names

def mask_sequence_name(sequence_name, env_from, env_to, mgnifams_fasta_dict):
    masked_name = f"{sequence_name}/{env_from}_{env_to}"

    if sequence_name in mgnifams_fasta_dict:
        sub_sequence = str(mgnifams_fasta_dict[sequence_name].seq[env_from - 1:env_to])  # -1 because Python uses 0-based indexing
        masked_seq_record = SeqRecord(Seq(sub_sequence), id=masked_name, description="")

        write_fasta_sequences(masked_seq_record, tmp_family_sequences_path, "a")

    return masked_name

def filter_recruited(evalue_threshold, length_threshold, mgnifams_fasta_dict, exit_flag):
    start_time = time.time()

    os.remove(tmp_family_sequences_path)
    filtered_sequences = []
    with open(tmp_domtblout_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split()
                evalue = float(columns[6])
                qlen = float(columns[5])
                env_from = int(columns[19])
                env_to = int(columns[20])
                env_length = env_to - env_from + 1
                if evalue < evalue_threshold:
                    if (exit_flag or env_length >= length_threshold * qlen): # only evalue filter when exiting, taking in shorter sequences
                        sequence_name = columns[0]
                        tlen = float(columns[2])
                        if (env_length < tlen):
                            sequence_name = mask_sequence_name(sequence_name, env_from, env_to, mgnifams_fasta_dict)
                        else:
                            write_fasta_sequences([mgnifams_fasta_dict[sequence_name]], tmp_family_sequences_path, "a")

                        filtered_sequences.append(sequence_name)
    
    with open(log_file, 'a') as file:
        file.write("filter_recruited: ")
        file.write(str(time.time() - start_time) + "\n")

    return filtered_sequences

def check_seed_membership(original_sequence_names, filtered_seq_names):
    def extract_first_part(sequence_name):
        return sequence_name.split('/')[0]

    original_first_parts = set(map(extract_first_part, original_sequence_names))
    filtered_first_parts = set(map(extract_first_part, filtered_seq_names))
    common_first_parts = original_first_parts & filtered_first_parts
    common_count = len(common_first_parts)
    percentage_membership = common_count / len(original_sequence_names)

    return percentage_membership

def run_hmmalign():
    start_time = time.time()

    hmmalign_command = ["hmmalign", "-o", tmp_align_msa_path, "--amino", tmp_hmm_path, tmp_family_sequences_path]
    subprocess.run(hmmalign_command, stdout=subprocess.DEVNULL)

    with open(log_file, 'a') as file:
        file.write("run_hmmalign: ")
        file.write(str(time.time() - start_time) + "\n")

def get_sequences_from_stockholm(file):
    with open(file, "r") as f:
        alignment = AlignIO.read(f, "stockholm")

    return alignment

def run_esl_weight(threshold=0.8):
    start_time = time.time()

    shutil.copy(tmp_align_msa_path, tmp_intermediate_esl_path)
    with open(log_file, 'a') as file:
        file.write("Files moved to tmp_intermediate_esl_path; ")

    while True:
        esl_weight_command = ["esl-weight", "--amino", "-f", "--idf", str(threshold), "-o", tmp_esl_weight_path, tmp_intermediate_esl_path]
        subprocess.run(esl_weight_command, stdout=subprocess.DEVNULL)
        with open(log_file, 'a') as file:
            file.write("esl_weight_command subprocess finished.\n")
        number_of_remaining_sequences = len(get_sequences_from_stockholm(tmp_esl_weight_path))
        with open(log_file, 'a') as file:
            file.write("Remaining sequences: " + str(number_of_remaining_sequences) + "\n")
        if (number_of_remaining_sequences <= 2000):
            break
        else:
            threshold -= 0.1
            with open(log_file, 'a') as file:
                file.write("Rerunning esl_weight; ")
            shutil.copy(tmp_esl_weight_path, tmp_intermediate_esl_path)

    with open(log_file, 'a') as file:
        file.write("run_esl_weight: ")
        file.write(str(time.time() - start_time) + "\n")

def extract_RF():
    with open(tmp_seed_msa_sto_path, 'r') as file:
        # Extract lines starting with "#=GC RF"
        relevant_lines = [line.strip() for line in file if line.startswith("#=GC RF")]

    # Keep only 'x's and '.'
    cleaned_lines = [''.join(filter(lambda c: c == 'x' or c == '.', line)) for line in relevant_lines]

    # Combine lines into a single sequence
    combined_sequence = ''.join(cleaned_lines)

    with open(tmp_rf_path, 'w') as output_file:
        output_file.write(combined_sequence)
        
def filter_out_redundant():
    # Step 1: Read identifiers from tmp_esl_weight_path using AlignIO
    identifiers = set()
    alignment = get_sequences_from_stockholm(tmp_esl_weight_path)
    for record in alignment:
        identifiers.add(record.id)

    # Step 2: Read and filter sequences in tmp_family_sequences_path
    filtered_sequences = []
    kept_identifiers = [] 
    for record in SeqIO.parse(tmp_family_sequences_path, "fasta"):
        if record.id in identifiers:
            filtered_sequences.append(record)
            kept_identifiers.append(record.id.split('/')[0])

    # Step 3: Write the filtered sequences back
    write_fasta_sequences(filtered_sequences, tmp_family_sequences_path, "w")

    return kept_identifiers

def append_family_file(iteration, family_members):
    lines = [f"{iteration}\t{member}\n" for member in family_members]
    with open(updated_refined_families_tsv_file, 'a') as file:
        file.writelines(lines)

    with open(log_file, 'a') as file:
        file.write(f"E: {iteration}, s: {len(family_members)}\n")

def parse_protein_region(protein_id):
    number_of_underscores = protein_id.count('_')
    if (number_of_underscores == 0):
        protein_rep = protein_id
        region      = "-"
    elif (number_of_underscores == 1):
        parts       = protein_id.split('/')
        protein_rep = parts[0]
        region      = parts[1].replace("_", "-")
    elif (number_of_underscores == 2):
        parts       = protein_id.split('_')
        protein_rep = parts[0]
        region      = f"{parts[1]}-{parts[2]}"
    elif (number_of_underscores == 3):
        parts        = protein_id.split('_')
        protein_rep  = parts[0]
        start        = int(parts[1])
        region_parts = protein_id.split('/')[1].split('_')
        region       = f"{start + int(region_parts[0]) - 1}-{start + int(region_parts[1]) - 1}"

    return protein_rep, region

def extract_first_stockholm_sequence():
    alignment = AlignIO.read(tmp_align_msa_path, "stockholm")
    first_record = alignment[0]
    protein_rep, region = parse_protein_region(first_record.id)
    return len(alignment), protein_rep, region 

def append_family_metadata(iteration):
    family_size, protein_rep, region = extract_first_stockholm_sequence()
    with open(updated_family_metadata_file, 'a') as file:
        file.writelines(f"{iteration},{family_size},{protein_rep},{region}\n")

def move_produced_models(iteration):
    shutil.move(tmp_seed_msa_sto_path, os.path.join(seed_msa_folder, f'{iteration}.sto'))
    shutil.move(tmp_align_msa_path, os.path.join(align_msa_folder, f'{iteration}.sto'))
    shutil.move(tmp_hmm_path, os.path.join(hmm_folder, f'{iteration}.hmm'))
    shutil.move(tmp_domtblout_path, os.path.join(domtblout_folder, f'{iteration}.domtblout'))
    shutil.move(tmp_rf_path, os.path.join(rf_folder, f'{iteration}.txt'))

def get_final_family_original_names(filtered_seq_names):
    family_members = {name.split('/')[0] for name in filtered_seq_names}
    return family_members

def get_cluster_reps(clusters_bookkeeping_df, family_members):
    # Match family_members to their family representatives and keep unique
    unique_reps = clusters_bookkeeping_df[clusters_bookkeeping_df['member'].isin(family_members)].index.unique()

    return unique_reps

def update_clusters_bookkeeping_df(clusters_bookkeeping_df, unique_reps):
    start_time = time.time()    

    # Flatten the 'member' column and filter based on unique_reps, then get unique values
    sequences_to_remove = clusters_bookkeeping_df.loc[clusters_bookkeeping_df.index.isin(unique_reps), 'member']
    if isinstance(sequences_to_remove.iloc[0], list):  # Check if the 'member' column contains lists
        sequences_to_remove = sequences_to_remove.explode().unique().tolist()
    else:
        sequences_to_remove = sequences_to_remove.unique().tolist()
    # Remove all lines having these family_reps from the clusters_bookkeeping_df
    indices_to_remove = clusters_bookkeeping_df.index[clusters_bookkeeping_df.index.isin(unique_reps)]
    clusters_bookkeeping_df.drop(index=indices_to_remove, inplace=True)

    with open(log_file, 'a') as file:
        file.write("update_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df, sequences_to_remove

def update_mgnifams_fasta_dict(mgnifams_fasta_dict, sequences_to_remove):
    start_time = time.time()
    
    # Convert sequences_to_remove to a set for faster lookup
    sequences_to_remove_set = set(sequences_to_remove)
    # Remove sequences in place
    for seq in sequences_to_remove_set:
        mgnifams_fasta_dict.pop(seq, None)

    write_fasta_sequences(mgnifams_fasta_dict.values(), updated_mgnifams_input_fasta_file, "w")

    with open(log_file, 'a') as file:
        file.write("update_mgnifams_fasta_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return mgnifams_fasta_dict

def remove_tmp_files():
    for item in os.listdir(tmp_folder):
        item_path = os.path.join(tmp_folder, item)
        if os.path.isfile(item_path):
            os.remove(item_path)

def main():
    parse_args()
    define_globals()

    copy_updated_input_files()
    clusters_bookkeeping_df = load_clusters_bookkeeping_df()
    clusters_bookkeeping_df = exclude_existing_families(clusters_bookkeeping_df) # restarting mechanism
    clusters_bookkeeping_df = exclude_discarded_clusters(clusters_bookkeeping_df) # restarting mechanism
    mgnifams_fasta_dict = create_mgnifams_fasta_dict()
    iteration = arg_iteration
    while True:
        iteration += 1
        largest_family_rep, family_members = get_next_largest_family(clusters_bookkeeping_df)
        if not family_members or len(family_members) < arg_minimum_family_members:
            with open(log_file, 'a') as file:
                file.write("Exiting all...")
            break
        
        write_family_fasta_file(family_members, mgnifams_fasta_dict)
        total_checked_sequences = family_members
        discard_flag = False
        exit_flag = False
        family_iteration = 0
        filtered_seq_names = []
        while True:
            family_iteration += 1
            if (family_iteration > 3):
                exit_flag = True

            with open(log_file, 'a') as file:
                if (exit_flag):
                    file.write("Exiting-3 loops.\n")
                file.write(str(family_iteration) + "\n")
            
            run_msa(len(family_members))
            original_sequence_names = trim_seed_msa()

            if not exit_flag: # main strategy branch
                run_hmmbuild(tmp_seed_msa_path, ["--amino", "--informat", "afa"])
                run_hmmsearch()
                recruited_sequence_names = extract_sequence_names_from_domtblout()
                filtered_seq_names = filter_recruited(evalue_threshold, length_threshold, mgnifams_fasta_dict, exit_flag) # also writes in tmp_family_sequences_path
                if (len(filtered_seq_names) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    break
                run_hmmalign()
                length_seqs_for_esl = len(get_sequences_from_stockholm(tmp_align_msa_path))
                if (length_seqs_for_esl > 70000):
                    discard_flag = True
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} too many sequences for esl ({length_seqs_for_esl}).\n")
                    break
                new_recruited_sequences = set(recruited_sequence_names) - set(total_checked_sequences)
                with open(log_file, 'a') as file:
                    file.write("new_recruited_sequences calculated.\n")
                if not new_recruited_sequences:
                    exit_flag = True
                    with open(log_file, 'a') as file:
                        file.write("Exiting-CONVERGED: no new sequences recruited.\n")
                    with open(updated_converged_families_file, 'a') as file:
                        file.write(f"{iteration}\n")

            if exit_flag: # exit strategy branch
                with open(log_file, 'a') as file:
                    file.write("Exiting branch strategy:\n")
                run_hmmbuild(tmp_seed_msa_path, ["-O", tmp_seed_msa_sto_path])
                run_hmmbuild(tmp_seed_msa_sto_path, ["--hand"])
                extract_RF()
                run_hmmsearch()
                filtered_seq_names = filter_recruited(evalue_threshold, length_threshold, mgnifams_fasta_dict, exit_flag) # also writes in tmp_family_sequences_path
                if (len(filtered_seq_names) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    break

                membership_percentage = check_seed_membership(original_sequence_names, filtered_seq_names)
                if (membership_percentage < 0.9):
                    discard_flag = True
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} seed percentage in MSA is {membership_percentage}\n")
                    break
                elif (membership_percentage < 1):
                    with open(log_file, 'a') as file:
                        file.write(f"Warning: {iteration} seed percentage in MSA is {membership_percentage}\n")
                
                run_hmmalign()
                break

            # main strategy branch continue
            total_checked_sequences += list(new_recruited_sequences)
            with open(log_file, 'a') as file:
                file.write("total_checked_sequences calculated and starting run_esl_weight\n")
            run_esl_weight()
            family_members = filter_out_redundant() # also writes in tmp_family_sequences_path

        # Exiting family loop
        if (discard_flag): # unsuccessfully
            with open(log_file, 'a') as file:
                file.write("Discarding cluster " + str(largest_family_rep) + "\n")

            iteration -= 1
            unique_reps = [largest_family_rep]
            with open(updated_discarded_clusters_file, 'a') as outfile:
                outfile.write(unique_reps[0] + "\n")
        else: # successfully
            append_family_file(iteration, filtered_seq_names)
            append_family_metadata(iteration)
            move_produced_models(iteration)
            # prepare sequences to remove
            family_original_names = get_final_family_original_names(filtered_seq_names)
            unique_reps = get_cluster_reps(clusters_bookkeeping_df, family_original_names)

        clusters_bookkeeping_df, sequences_to_remove = update_clusters_bookkeeping_df(clusters_bookkeeping_df, unique_reps)
        mgnifams_fasta_dict = update_mgnifams_fasta_dict(mgnifams_fasta_dict, sequences_to_remove)
        remove_tmp_files()

    # End of all families

if __name__ == "__main__":
    main()
