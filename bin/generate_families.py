#!/usr/bin/env python3

import argparse
import os
import time
import re
import shutil

import pandas as pd
import pyfastx
import pyfamsa
import pyhmmer
import pytrimal
import numpy as np

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
        converged_families_file, family_metadata_file, family_reps_file, \
        tmp_folder, seed_msa_folder, \
        align_msa_folder, hmm_folder, \
        domtblout_folder, rf_folder, \
        evalue_threshold, length_threshold, \
        tmp_family_sequences_path, tmp_seed_msa_path, tmp_align_msa_path, \
        tmp_hmm_path, tmp_domtblout_path, \
        tmp_sequences_to_remove_path, tmp_rf_path

    logs_folder                = "logs"
    refined_families_folder    = "refined_families"
    discarded_clusters_folder  = "discarded_clusters"
    successful_clusters_folder = "successful_clusters"
    converged_families_folder  = "converged_families"
    family_metadata_folder     = "family_metadata"
    family_reps_folder         = "family_reps"
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
        successful_clusters_folder, converged_families_folder, family_metadata_folder, family_reps_folder
    ]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    tmp_family_sequences_path    = os.path.join(tmp_folder, 'family_sequences.fa')
    tmp_seed_msa_path            = os.path.join(tmp_folder, 'seed_msa.sto')
    tmp_align_msa_path           = os.path.join(tmp_folder, 'align_msa.sto')
    tmp_hmm_path                 = os.path.join(tmp_folder, 'model.hmm')
    tmp_domtblout_path           = os.path.join(tmp_folder, 'domtblout.txt')
    tmp_sequences_to_remove_path = os.path.join(tmp_folder, 'sequences_to_remove.txt')
    tmp_rf_path                  = os.path.join(tmp_folder, 'rf.txt')

    log_file                  = os.path.join(logs_folder               , f'{arg_chunk_num}.txt')
    refined_families_tsv_file = os.path.join(refined_families_folder   , f'{arg_chunk_num}.tsv')
    discarded_clusters_file   = os.path.join(discarded_clusters_folder , f'{arg_chunk_num}.csv')
    successful_clusters_file  = os.path.join(successful_clusters_folder, f'{arg_chunk_num}.txt')
    converged_families_file   = os.path.join(converged_families_folder , f'{arg_chunk_num}.txt')
    family_metadata_file      = os.path.join(family_metadata_folder    , f'{arg_chunk_num}.csv')
    family_reps_file          = os.path.join(family_reps_folder        , f'{arg_chunk_num}.fasta')

def create_empty_output_files():
    open(refined_families_tsv_file, 'w').close()
    open(discarded_clusters_file  , 'w').close()
    open(successful_clusters_file , 'w').close()
    open(converged_families_file  , 'w').close()
    open(family_metadata_file     , 'w').close()
    open(family_reps_file         , 'w').close()

def load_clusters_df():
    return pd.read_csv(arg_clusters_chunk, sep='\t', header=None, names=['representative', 'member'], dtype=str)

def log_time(start_time, text, mode='a'):
    with open(log_file, mode) as file:
        file.write(text)
        file.write(str(time.time() - start_time) + "\n")

def create_mgnifams_pyfastx_obj():
    start_time = time.time()

    mgnifams_pyfastx_obj = pyfastx.Fasta(arg_mgnifams_input_fasta_file)

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

def get_fasta_sequences(pyfastx_obj, headers):
    """Retrieve sequences from pyfastx_obj given a list of headers."""
    return [(header, str(pyfastx_obj[header])) for header in headers if header in pyfastx_obj]

def write_fasta_sequences(sequences, file_path, mode):
    with open(file_path, mode) as output_handle:
        for header, sequence in sequences:
            output_handle.write(f">{header}\n{sequence}\n")

def write_family_fasta_file(members, pyfastx_obj):
    start_time = time.time()

    sequences_to_write = get_fasta_sequences(pyfastx_obj, members)
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

    log_time(start_time, "run_initial_msa (pyfamsa): ")

def run_hmmbuild(iteration, hand=False):
    """Runs HMMER's hmmbuild using pyhmmer."""
    start_time = time.time()                    

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.MSAFile(tmp_seed_msa_path, digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
    msa.name = f"{arg_chunk_num}_{iteration}".encode()

    architecture = "hand" if hand else "fast"
    builder = pyhmmer.plan7.Builder(alphabet, architecture=architecture)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)

    if (architecture == "hand"): # only write to outfile if last iteration, else using the built object
        with open(tmp_hmm_path, "wb") as output_file:
            hmm.write(output_file)

    log_time(start_time, "run_hmmbuild (pyhmmer): ")

    return hmm, hmm.consensus

def mask_sequence(sequence, env_from, env_to):
    # Unpack the tuple (header, seq)
    header, seq = sequence[0]  # sequence[0] is a tuple

    # Modify header and sequence
    header = f"{header}/{env_from}_{env_to}"
    seq = seq[env_from - 1:env_to]  # -1 for 0-based indexing

    # Replace the old tuple with the new one
    sequence[0] = (header, seq)

    return sequence

def run_hmmsearch(hmm, pyhmmer_seqs, pyfastx_obj, exit_flag): # hmm obj changes, pyhmmer_seqs is the same, whole protein library
    """Runs HMMER's hmmsearch using pyhmmer."""
    start_time = time.time()

    filtered_sequences = []
    os.remove(tmp_family_sequences_path) # emptying file

    for top_hits in pyhmmer.hmmer.hmmsearch(hmm, pyhmmer_seqs, cpus=int(arg_cpus), E=evalue_threshold):
        with open(tmp_domtblout_path, "wb") as fh:
            top_hits.write(fh, format="domains") # --domtblout
        qlen = top_hits.query.M
        for hit in top_hits:
            sequence_name = hit.name.decode()
            # print(f"Target name: {hit.name.decode()}, Hit E-value: {hit.evalue}") # debug
            tlen = hit.length
            for subhit in hit.domains:
                # print(f"TARGET SEQUENCE: {subhit.alignment.target_sequence}") # debug, different than target_env_sequence that we want, not provided by pyhmmer
                env_length = subhit.env_to - subhit.env_from + 1
                # print(f"env from-to: {subhit.env_from}-{subhit.env_to}") # debug
                if (exit_flag or env_length >= length_threshold * qlen):
                    sequence = get_fasta_sequences(pyfastx_obj, [sequence_name])
                    if (env_length < tlen):
                        sequence = mask_sequence(sequence, subhit.env_from, subhit.env_to)
                    write_fasta_sequences(sequence, tmp_family_sequences_path, "a")
                    filtered_sequences.append(sequence[0][0]) # [0][0] is the header (list of tuples)
                    # print(f"ENV SEQUENCE: {sequence[0][1]}") # debug

    log_time(start_time, "run_hmmsearch (pyhmmer): ")

    return filtered_sequences

def run_hmmalign(hmm):
    """Runs HMMER's hmmalign using pyhmmer."""
    start_time = time.time()

    num_seqs_result = 0
    non_gap_seq_length = 0

    with pyhmmer.easel.SequenceFile(tmp_family_sequences_path, digital=True) as seq_file:
        seqs = seq_file.read_block()
        hmmalign_res = pyhmmer.hmmer.hmmalign(hmm, seqs, trim=True)

        with open(tmp_align_msa_path, "wb") as outfile: # full MSA written out here
            hmmalign_res.write(outfile, format="stockholm") # expected ['stockholm', 'pfam', 'a2m', 'psiblast', 'selex', 'afa', 'clustal', 'clustallike', 'phylip' or 'phylips']

        num_seqs_result = len(hmmalign_res.names)
        non_gap_seq_length = len(re.sub(r"[.-]", "", hmmalign_res.alignment[0]))

    log_time(start_time, "run_hmmalign (pyhmmer): ")

    return num_seqs_result, non_gap_seq_length

def unmask_sequence_names(names):
    return [name.split('/')[0] for name in names]

def run_pytrimal_reps(threshold=0.8):
    start_time = time.time()

    # Load the MSA
    bio_ali = AlignIO.read(tmp_align_msa_path, "stockholm")
    names = [record.id.encode() for record in bio_ali]
    sequence = [bytes(record.seq) for record in bio_ali]
    msa = pytrimal.Alignment(names, sequence)

    # and remove redundant sequences
    while True:
        # print(len(list(msa.names))) # debug
        repTrimmer = pytrimal.RepresentativeTrimmer(identity_threshold=threshold)
        msa = repTrimmer.trim(msa)

        with open(log_file, 'a') as file:
            file.write("run_pytrimal_reps finished.\n")
        number_of_remaining_sequences = len(list(msa.names))
        with open(log_file, 'a') as file:
            file.write("Remaining sequences: " + str(number_of_remaining_sequences) + "\n")

        if (number_of_remaining_sequences <= 2000): # TODO pytrimal -clusters 2000 max once, if/when fixed
            break
        else:
            threshold -= 0.1
            with open(log_file, 'a') as file:
                file.write("Rerunning run_pytrimal_reps; ")

    log_time(start_time, "run_pytrimal_reps: ")

    return msa

def write_filtered_sto_to_seed_msa_file(name_set, start_pos, end_pos):
    with open(tmp_align_msa_path, 'r') as infile, open(tmp_seed_msa_path, 'w') as outfile:
        for line in infile:
            # Split line by spaces
            split_line = line.split()

            # Check if the line meets any of the specified conditions
            if not split_line or len(split_line) == 1:  # Line is empty or //
                outfile.write(line)
            elif split_line[1] == 'STOCKHOLM':  # The second split is 'STOCKHOLM'
                outfile.write(line)
            elif split_line[0] in name_set:  # First split is in name_set
                original_substring = split_line[1]  # The part to modify
                modified_substring = original_substring[start_pos:end_pos]  # Extract substring
                line = line.replace(original_substring, modified_substring, 1)
                outfile.write(line)
            elif split_line[1] in name_set:  # Second split is in name_set
                original_substring = split_line[3]  # The part to modify
                modified_substring = original_substring[start_pos:end_pos]  # Extract substring
                line = line.replace(original_substring, modified_substring, 1)
                outfile.write(line)
            elif split_line[0] == '#=GC':  # First split is '#=GC'
                original_substring = split_line[2]  # The part to modify
                modified_substring = original_substring[start_pos:end_pos]  # Extract substring
                line = line.replace(original_substring, modified_substring, 1)
                outfile.write(line)
                if split_line[1] == 'RF':
                    start_pos = 0
                    end_pos = end_pos - 200 # multi-line sto MSA, 200 aa per line

def calculate_trim_positions(sequence_matrix, occupancy_threshold):
    numeric_matrix = np.where(sequence_matrix == '-', 0, 1)
    num_rows = numeric_matrix.shape[0]
    column_sums = np.sum(numeric_matrix, axis=0)
    column_sums_percentage = column_sums / num_rows
    start_position = np.argmax(column_sums_percentage > occupancy_threshold)
    end_position = len(column_sums_percentage) - np.argmax(column_sums_percentage[::-1] > occupancy_threshold) - 1

    return start_position, end_position

def clip_ends(msa, occupancy_threshold=0.5):
    sequence_matrix = np.array([list(seq) for seq in msa.sequences])
    start_position, end_position = calculate_trim_positions(sequence_matrix, occupancy_threshold)
    
    name_set = {name.decode() for name in msa.names}
    write_filtered_sto_to_seed_msa_file(name_set, start_position, end_position)

# def get_anchor_points(main_seq, sub_seq):
#     match = re.search(sub_seq, main_seq) # TODO be careful if there is ever a multiple match, will get wrong anchor points for the MSA (try re.finditer)
#     start_position, end_position = match.span()

#     return start_position, end_position

# def run_pytrimal_terminalonly(msa): # TODO switch to this when and if terminal_only bug fixed
#     trimmer = pytrimal.ManualTrimmer(gap_threshold=0.5)
#     trimmed_alignment = trimmer.trim(msa).terminal_only() # trim ends only

#     start_position, end_position = get_anchor_points(msa.sequences[0], trimmed_alignment.sequences[0])

#     name_set = {name.decode() for name in trimmed_alignment.names}
#     write_filtered_sto_to_seed_msa_file(name_set, start_position, end_position)

#     # TODO manipulate the object below, Bateman trim?
#     # for name, aligned in zip(hmmalign_res.names, hmmalign_res.alignment):
#     #     print(name.decode(), " ", aligned)

def extract_RF():
    with open(tmp_seed_msa_path, 'r') as file:
        # Extract lines starting with "#=GC RF"
        relevant_lines = [line.strip() for line in file if line.startswith("#=GC RF")]

    # Keep only 'x's and '.'
    cleaned_lines = [''.join(filter(lambda c: c == 'x' or c == '.', line)) for line in relevant_lines]

    # Combine lines into a single sequence
    combined_sequence = ''.join(cleaned_lines)

    with open(tmp_rf_path, 'w') as output_file:
        output_file.write(combined_sequence)

def check_seed_membership(original_sequence_names, filtered_seq_names):
    def extract_first_part(sequence_name):
        return sequence_name.split('/')[0]

    original_first_parts = set(map(extract_first_part, original_sequence_names))
    filtered_first_parts = set(map(extract_first_part, filtered_seq_names))
    common_first_parts = original_first_parts & filtered_first_parts
    common_count = len(common_first_parts)
    percentage_membership = common_count / len(original_sequence_names)

    return percentage_membership

def parse_protein_region(protein_id): # TODO check if works, update or remove
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

def parse_protein_name(seq_name, seq_length, seq_whole_name, original_length, start, end):

    if ( (end-start) == original_length):
        return seq_name
    else:
        splits = seq_name.split('_')
        old_length = len(seq_whole_name)
        if (len(splits) == 3):
            new_start = start + int(splits[1])
            new_end = new_start + seq_length -1
            new_name = splits[0] + "/" + str(new_start) + "-" + str(new_end)
        else:
            new_name = splits[0] + "/" + str(start + 1) + "-" + str(end)

        new_name = new_name.ljust(old_length)  # Pad with spaces to match original length

        return new_name

def renumber_seed_sto_msa(in_sto_file, out_sto_file, pyfastx_obj):
    with open(in_sto_file, 'r') as infile, open(out_sto_file, 'w') as outfile:

        previous_seq_name = ""
        for line in infile:
            # Split line by spaces
            split_line = line.split()

            # Check if the line meets any of the specified conditions
            if not split_line or len(split_line) == 1 or split_line[0] == '#=GC':  # Line is empty or // or '#=GC'
                outfile.write(line)
            elif split_line[1] == 'STOCKHOLM':  # The second split is 'STOCKHOLM'
                outfile.write(line)
            elif split_line[0] != '#=GR' and split_line[0] != '#=GC':  # First split is first protein encounter
                seq_name = split_line[0].split("/")[0]
                seq = re.sub(r"[.-]", "", split_line[1]).upper()
                original_seq = get_fasta_sequences(pyfastx_obj, [seq_name])[0][1]
                start = original_seq.find(seq)
                end = start + len(seq)
                previous_seq_name = parse_protein_name(seq_name, len(seq), split_line[0], len(original_seq), start, end)
                line = line.replace(split_line[0], previous_seq_name, 1)
                outfile.write(line)
            elif split_line[0] == '#=GR':
                line = line.replace(split_line[1], previous_seq_name, 1)
                outfile.write(line)

def renumber_full_sto_msa_and_write_tsv_metadata(in_sto_file, pyfastx_obj, arg_chunk_num, iteration, \
                                                full_msa_num_seqs, consensus):
    
    out_sto_file = os.path.join(align_msa_folder, f'{arg_chunk_num}_{iteration}.sto')
    with open(in_sto_file, 'r') as infile, open(refined_families_tsv_file, 'a') as familyfile, \
        open(family_metadata_file, 'a') as metadatafile, open(family_reps_file, 'a') as repsfile, \
        open(out_sto_file, 'w') as outfile:

        previous_seq_name = ""
        rep_flag = True
        for line in infile:
            # Split line by spaces
            split_line = line.split()

            # Check if the line meets any of the specified conditions
            if not split_line or len(split_line) == 1 or split_line[0] == '#=GC':  # Line is empty or // or '#=GC'
                outfile.write(line)
            elif split_line[1] == 'STOCKHOLM':  # The second split is 'STOCKHOLM'
                outfile.write(line)
            elif split_line[0] != '#=GR' and split_line[0] != '#=GC':  # First split is first protein encounter
                seq_name = split_line[0].split("/")[0]
                seq = re.sub(r"[.-]", "", split_line[1]).upper()
                original_seq = get_fasta_sequences(pyfastx_obj, [seq_name])[0][1]
                start = original_seq.find(seq)
                end = start + len(seq)
                previous_seq_name = parse_protein_name(seq_name, len(seq), split_line[0], len(original_seq), start, end)
                line = line.replace(split_line[0], previous_seq_name, 1)
                outfile.write(line)
                familyfile.write(f"{iteration}\t{previous_seq_name}\n")
                if (rep_flag):
                    metadatafile.write(f"{iteration},{full_msa_num_seqs},\"{previous_seq_name.strip()}\",{len(seq)},{seq},{consensus}\n")
                    repsfile.write(f">{previous_seq_name.strip()}\t{arg_chunk_num}_{iteration}\n{seq}\n")
                    rep_flag = False
            elif split_line[0] == '#=GR':
                line = line.replace(split_line[1], previous_seq_name, 1)
                outfile.write(line)

def append_family_metadata(protein_rep, region, sequence, iteration, full_msa_num_seqs, consensus):  # TODO maybe along with full MSA
    with open(family_metadata_file, 'a') as file:
        file.writelines(f"{iteration},{full_msa_num_seqs},\"{protein_rep}\",{region},{len(sequence)},{sequence},{consensus}\n")
    with open(family_reps_file, 'a') as file:
        file.writelines(
        f">{protein_rep}\t{arg_chunk_num}_{iteration}\n{sequence}\n"
        if region == "-" 
        else f">{protein_rep}/{region}\t{arg_chunk_num}_{iteration}\n{sequence}\n"
    )

def move_produced_models(iteration):
    shutil.move(tmp_hmm_path,      os.path.join(hmm_folder,       f'{arg_chunk_num}_{iteration}.hmm'))
    shutil.move(tmp_domtblout_path,os.path.join(domtblout_folder, f'{arg_chunk_num}_{iteration}.domtblout'))
    shutil.move(tmp_rf_path,       os.path.join(rf_folder,        f'{arg_chunk_num}_{iteration}.txt'))

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
        
        original_sequence_names = family_members
        write_family_fasta_file(family_members, mgnifams_pyfastx_obj)

        run_initial_msa(family_members, mgnifams_pyfastx_obj)

        total_checked_sequences = []
        filtered_seq_names = []
        full_msa_num_seqs = 0
        discard_flag = False
        discard_reason = ""
        discard_value = 0.0
        exit_flag = False
        hmm = ""
        consensus = ""
        family_iteration = 0
        while True:
            family_iteration += 1
            if (family_iteration > 3):
                exit_flag = True

            with open(log_file, 'a') as file:
                if (exit_flag):
                    file.write("Exiting-3 loops.\n")
                file.write(str(family_iteration) + "\n")

            if not exit_flag: # main strategy branch
                hmm, _ = run_hmmbuild(iteration)

                filtered_seq_names = run_hmmsearch(hmm, pyhmmer_seqs, mgnifams_pyfastx_obj, exit_flag)
                if (len(filtered_seq_names) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    discard_reason = "low complexity model - confounding cluster"
                    discard_value = 0.0
                    break

                new_recruited_sequences = set(unmask_sequence_names(filtered_seq_names)) - set(total_checked_sequences) # new_recruited_sequences always has something at first turn since total_checked_sequences starts empty []
                total_checked_sequences += list(new_recruited_sequences)

                if not new_recruited_sequences:
                    exit_flag = True
                    with open(log_file, 'a') as file:
                        file.write("Exiting-CONVERGED: no new sequences recruited.\n")
                    with open(converged_families_file, 'a') as file:
                        file.write(f"{iteration}\n")

            if exit_flag: # exit strategy branch
                with open(log_file, 'a') as file:
                    file.write("Exiting branch strategy:\n")

                _, consensus = run_hmmbuild(iteration, hand=True)
                extract_RF()

                filtered_seq_names = run_hmmsearch(hmm, pyhmmer_seqs, mgnifams_pyfastx_obj, exit_flag)
                if (len(filtered_seq_names) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    discard_reason = "low complexity model - confounding cluster"
                    discard_value = 0.0
                    break

                membership_percentage = check_seed_membership(original_sequence_names, unmask_sequence_names(filtered_seq_names))
                if (membership_percentage < 0.9):
                    discard_flag = True
                    discard_reason = "few seed sequences remained"
                    discard_value = membership_percentage
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} seed percentage in MSA is {membership_percentage}\n")
                    break
                elif (membership_percentage < 1):
                    with open(log_file, 'a') as file:
                        file.write(f"Warning: {iteration} seed percentage in MSA is {membership_percentage}\n")
                
                full_msa_num_seqs, non_gap_seq_length = run_hmmalign(hmm) # final full MSA, including smaller sequences
                if (non_gap_seq_length < 100):
                    discard_flag = True
                    discard_reason = "family representative length too small"
                    discard_value = non_gap_seq_length
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} rep length is only {non_gap_seq_length}\n")
                    break
                elif (non_gap_seq_length > 2000):
                    discard_flag = True
                    discard_reason = "family representative length too large"
                    discard_value = non_gap_seq_length
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} rep length is over {non_gap_seq_length}\n")
                    break

                break # break from main strategy

            # main strategy continue, if not converged
            run_hmmalign(hmm)
            msa = run_pytrimal_reps() # removes redundant sequences
            clip_ends(msa) # removes gaps above threshold at ends
            # run_pytrimal_terminalonly(msa) # removes gaps above threshold at ends

        # Exiting family loop
        if (discard_flag): # unsuccessfully
            with open(log_file, 'a') as file:
                file.write("Discarding cluster " + family_rep + "\n")
            
            with open(discarded_clusters_file, 'a') as outfile:
                outfile.write(str(family_rep) + "," + discard_reason + "," + str(discard_value) + "\n")
            iteration -= 1 # keep proper track of family ids
        else: # successfully
            with open(successful_clusters_file, 'a') as outfile:
                outfile.write(str(family_rep) + "\n")
            
            renumber_seed_sto_msa(tmp_seed_msa_path, os.path.join(seed_msa_folder, f'{arg_chunk_num}_{iteration}.sto'), mgnifams_pyfastx_obj)
            renumber_full_sto_msa_and_write_tsv_metadata(tmp_align_msa_path, mgnifams_pyfastx_obj, \
                                                        arg_chunk_num, iteration, full_msa_num_seqs, consensus)
            move_produced_models(iteration)

        remove_tmp_files()

    # # End of all families
    with open(log_file, 'a') as file:
        file.write("DONE.")

if __name__ == "__main__":
    main()
