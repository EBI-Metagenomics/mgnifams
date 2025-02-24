#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import time
import pyfastx
import pyfamsa
import pyhmmer

import subprocess
import shutil
import numpy as np
from Bio.Seq import Seq # TODO remove
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO


def parse_args(args=None):
    global arg_clusters_chunk, arg_mgnifams_input_fasta_file, arg_cpus, arg_chunk_num
    parser = argparse.ArgumentParser(description="Process clustering data and extract sequences.")
    parser.add_argument("-c", "--clusters_chunk", required=True, type=str, help="Path to the clusters chunk file.")
    parser.add_argument("-f", "--fasta_file", required=True, type=str, help="Path to the MGnifams input FASTA file.")
    parser.add_argument("-p", "--cpus", required=True, type=str, help="Number of CPUs to use.")
    parser.add_argument("-n", "--chunk_num", required=True, type=str, help="Chunk number used for naming output files.")
    args = parser.parse_args(args)
    
    arg_clusters_chunk = args.clusters_chunk
    arg_mgnifams_input_fasta_file = args.fasta_file
    arg_cpus = args.cpus
    arg_chunk_num = args.chunk_num

def define_globals():
    global log_file, refined_families_tsv_file, \
        discarded_clusters_file, successful_clusters_file, \
        converged_families_file, family_metadata_file, \
        tmp_folder, seed_msa_folder, \
        align_msa_folder, hmm_folder, \
        domtblout_folder, rf_folder, \
        evalue_threshold, length_threshold, \
        tmp_family_sequences_path, tmp_seed_msa_path, \
        tmp_seed_msa_sto_path, tmp_align_msa_path, \
        tmp_hmm_path, tmp_domtblout_path, \
        tmp_intermediate_esl_path, tmp_esl_weight_path, \
        tmp_sequences_to_remove_path, tmp_rf_path

    logs_folder                = "logs"
    refined_families_folder    = "refined_families"
    discarded_clusters_folder  = "discarded_clusters"
    successful_clusters_folder = "successful_clusters"
    converged_families_folder  = "converged_families"
    family_metadata_folder     = "family_metadata"
    tmp_folder                 = "tmp"
    seed_msa_folder            = "seed_msa_sto"
    align_msa_folder           = "msa_sto"
    hmm_folder                 = "hmm"
    domtblout_folder           = "domtblout"
    rf_folder                  = "rf"
    evalue_threshold           = 0.001
    length_threshold           = 0.8

    for folder in [tmp_folder, seed_msa_folder, align_msa_folder, hmm_folder, domtblout_folder, rf_folder, \
        logs_folder, refined_families_folder, discarded_clusters_folder, \
        successful_clusters_folder, converged_families_folder, family_metadata_folder
    ]:
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

    log_file                  = os.path.join(logs_folder               , f'{arg_chunk_num}.txt')
    refined_families_tsv_file = os.path.join(refined_families_folder   , f'{arg_chunk_num}.tsv')
    discarded_clusters_file   = os.path.join(discarded_clusters_folder , f'{arg_chunk_num}.txt')
    successful_clusters_file  = os.path.join(successful_clusters_folder, f'{arg_chunk_num}.txt')
    converged_families_file   = os.path.join(converged_families_folder , f'{arg_chunk_num}.txt')
    family_metadata_file      = os.path.join(family_metadata_folder    , f'{arg_chunk_num}.csv')

def create_empty_output_files():
    open(refined_families_tsv_file, 'w').close()
    open(discarded_clusters_file  , 'w').close()
    open(successful_clusters_file , 'w').close()
    open(converged_families_file  , 'w').close()
    open(family_metadata_file     , 'w').close()

def load_clusters_df():
    return pd.read_csv(arg_clusters_chunk, sep='\t', header=None, names=['representative', 'member'], dtype=str)

def log_time(start_time, text, mode='a'):
    with open(log_file, mode) as file:
        file.write(text)
        file.write(str(time.time() - start_time) + "\n")

def create_mgnifams_pyfastx_obj():
    start_time = time.time()

    mgnifams_pyfastx_obj = pyfastx.Fasta(arg_mgnifams_input_fasta_file, key_func=str)

    log_time(start_time, "create_mgnifams_fasta_dict: ", 'w')

    return mgnifams_pyfastx_obj

def get_next_family(clusters_df):
    start_time = time.time()

    if clusters_df.empty:
        return None, None

    # Select the next family (first row in the dataFrame)
    next_family_rep = clusters_df.iloc[0]['representative']
    next_family_members = clusters_df.loc[clusters_df['representative'] == next_family_rep, 'member'].tolist()

    # Remove the selected family from the dataframe
    clusters_df.drop(clusters_df[clusters_df['representative'] == next_family_rep].index, inplace=True)

    # Logging execution time and results
    with open(log_file, 'a') as file:
        file.write(f"\nget_next_family: {time.time() - start_time}\n")
        file.write(f"S: {next_family_rep}, s: {len(next_family_members)}\n")
    
    return next_family_rep, next_family_members

def read_pyhmmer_seqs():
    start_time = time.time()

    # pre-fetching targets - fast, but needs the whole target database in memory
    with pyhmmer.easel.SequenceFile(arg_mgnifams_input_fasta_file, digital=True) as seq_file:
        seqs = seq_file.read_block()

    log_time(start_time, "read_pyhmmer_seqs (pyhmmer): ")

    return seqs

def write_fasta_sequences(sequences, file_path, mode):
    with open(file_path, mode) as output_handle:
        for header, sequence in sequences:
            output_handle.write(f">{header}\n{sequence}\n")

def write_family_fasta_file(members, pyfastx_obj):
    start_time = time.time()

    sequences_to_write = [(member, str(pyfastx_obj[member])) for member in members if member in pyfastx_obj]
    write_fasta_sequences(sequences_to_write, tmp_family_sequences_path, "w")

    log_time(start_time, "write_family_fasta_file (pyfastx): ")

def run_initial_msa(members, pyfastx_obj):
    """Aligns sequences using pyfamsa and writes the result to a FASTA file."""
    start_time = time.time()
    
    # Convert sequences to pyfamsa.Sequence objects
    sequences = [pyfamsa.Sequence(member.encode(), str(pyfastx_obj[member]).encode()) for member in members if member in pyfastx_obj]
    
    if not sequences:
        raise ValueError("No valid sequences found for alignment.")

    # Create an Aligner object
    aligner = pyfamsa.Aligner()
    
    # Perform the multiple sequence alignment
    alignment = aligner.align(sequences)

    # Write the aligned sequences to a FASTA file
    with open(tmp_seed_msa_path, "w") as output_handle:
        for gapped_seq in alignment:
            output_handle.write(f">{gapped_seq.id.decode()}\n{gapped_seq.sequence.decode()}\n")

    shutil.copy(tmp_seed_msa_path, tmp_align_msa_path) # copy to full alignment, in case hmmsearch doesn't recruit anything new at the first loop (never happens)

    log_time(start_time, "run_initial_msa (pyfamsa): ")

def run_hmmbuild(msa_file, chunk_num):
    """Runs HMMER's hmmbuild using pyhmmer."""
    start_time = time.time()                    

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.MSAFile(msa_file, digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
    msa.name = f"protein_family_{chunk_num}".encode()

    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)
    print(hmm.consensus) # TODO write to file and emit as output

    with open(tmp_hmm_path, "wb") as output_file:
        hmm.write(output_file)

    log_time(start_time, "run_hmmbuild (pyhmmer): ")

def run_hmmsearch(pyhmmer_seqs):
    """Runs HMMER's hmmsearch using pyhmmer."""
    start_time = time.time()

    with pyhmmer.plan7.HMMFile(tmp_hmm_path) as hmm:
        for top_hits in pyhmmer.hmmer.hmmsearch(hmm, pyhmmer_seqs, cpus=int(arg_cpus), E=evalue_threshold):
            print(len(top_hits))
            for hit in top_hits:
                print(f"Target: {hit.name.decode()}, E-value: {hit.evalue}")


    # Write the results to the output file
    # with open(tmp_domtblout_path, "w") as output_file:
    #     hits.write(output_file, format="domtbl")

    # hmmsearch_command = ["hmmsearch", "--cpu", arg_cpus, "--domtblout", tmp_domtblout_path, tmp_hmm_path, arg_mgnifams_input_fasta_file]
    # subprocess.run(hmmsearch_command, stdout=subprocess.DEVNULL)

    log_time(start_time, "run_hmmsearch (pyhmmer): ")

#=========================== UPDATED UP TO HERE ===========================================#

def run_msa(family_size):
    start_time = time.time()
    
    if (family_size > 1):
        mafft_command = ["mafft", "--quiet", "--retree", "2", "--maxiterate", "2", "--thread", "-1", tmp_family_sequences_path]
    else:
        with open(log_file, 'a') as file:
            file.write("Running 1-sequence version of mafft.\n")
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

def run_hmmbuild_hmmer(msa_file, extra_args): # TODO remove
    start_time = time.time()                    

    hmmbuild_command = ["hmmbuild", "--cpu", arg_cpus] + extra_args + [tmp_hmm_path, msa_file]
    subprocess.run(hmmbuild_command, stdout=subprocess.DEVNULL)
    
    with open(log_file, 'a') as file:
        file.write("run_hmmbuild: ")
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
    with open(refined_families_tsv_file, 'a') as file:
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
    with open(family_metadata_file, 'a') as file:
        file.writelines(f"{iteration},{family_size},{protein_rep},{region}\n")

def move_produced_models(iteration):
    shutil.move(tmp_seed_msa_sto_path, os.path.join(seed_msa_folder,  f'{arg_chunk_num}_{iteration}.sto'))
    shutil.move(tmp_align_msa_path,    os.path.join(align_msa_folder, f'{arg_chunk_num}_{iteration}.sto'))
    shutil.move(tmp_hmm_path,          os.path.join(hmm_folder,       f'{arg_chunk_num}_{iteration}.hmm'))
    shutil.move(tmp_domtblout_path,    os.path.join(domtblout_folder, f'{arg_chunk_num}_{iteration}.domtblout'))
    shutil.move(tmp_rf_path,           os.path.join(rf_folder,        f'{arg_chunk_num}_{iteration}.txt'))

def get_final_family_original_names(filtered_seq_names):
    family_members = {name.split('/')[0] for name in filtered_seq_names}
    return family_members

def remove_tmp_files():
    for item in os.listdir(tmp_folder):
        item_path = os.path.join(tmp_folder, item)
        if os.path.isfile(item_path):
            os.remove(item_path)

def main():
    parse_args()
    define_globals()
    create_empty_output_files()
    clusters_df = load_clusters_df()

    mgnifams_pyfastx_obj = create_mgnifams_pyfastx_obj()
    pyhmmer_seqs = read_pyhmmer_seqs()

    iteration = 0
    while True:
        iteration += 1
        family_rep, family_members = get_next_family(clusters_df)
        if not family_members:
            with open(log_file, 'a') as file:
                file.write("Exiting all...")
            break

        write_family_fasta_file(family_members, mgnifams_pyfastx_obj)

        run_initial_msa(family_members, mgnifams_pyfastx_obj)
        # TODO run trim_seed_msa / clipkit? Probably not needed
        # TODO renumber regions after clipping

        total_checked_sequences = family_members
        filtered_seq_names = []
        discard_flag = False
        exit_flag = False
        family_iteration = 0
        while True:
            family_iteration += 1
            if (family_iteration > 3):
                exit_flag = True

            with open(log_file, 'a') as file:
                if (exit_flag):
                    file.write("Exiting-3 loops.\n")
                file.write(str(family_iteration) + "\n")
            
            # run_msa(len(family_members)) # TODO remove from here?
            # original_sequence_names = trim_seed_msa() # TODO remove from here?

            if not exit_flag: # main strategy branch
                run_hmmbuild(tmp_seed_msa_path, arg_chunk_num)
                run_hmmsearch(pyhmmer_seqs)
                exit()

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
                    with open(converged_families_file, 'a') as file:
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
                file.write("Discarding cluster " + family_rep + "\n")
            
            with open(discarded_clusters_file, 'a') as outfile:
                outfile.write(str(family_rep) + "\n")
            iteration -= 1
        else: # successfully
            with open(successful_clusters_file, 'a') as outfile:
                outfile.write(str(family_rep) + "\n")
            append_family_file(iteration, filtered_seq_names)
            append_family_metadata(iteration)
            move_produced_models(iteration)
        
        remove_tmp_files()

    # # End of all families
    with open(log_file, 'a') as file:
        file.write("DONE.")

if __name__ == "__main__":
    main()
