#!/usr/bin/env python3

import argparse
import os
import time
import shutil
import subprocess

import pandas as pd
import pyfastx
import pyfamsa
import pyhmmer

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

    shutil.copy(tmp_seed_msa_path, tmp_align_msa_path) # copy to full alignment, in case hmmsearch doesn't recruit anything new at the first loop (never happens)

    log_time(start_time, "run_initial_msa (pyfamsa): ")

def run_hmmbuild(msa_file, hand=False): # TODO remove msa_file arg if always called for same global file
    """Runs HMMER's hmmbuild using pyhmmer."""
    start_time = time.time()                    

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.MSAFile(msa_file, digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
    msa.name = f"protein_family_{arg_chunk_num}".encode()

    architecture = "hand" if hand else "fast"
    builder = pyhmmer.plan7.Builder(alphabet, architecture=architecture)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)
    print(hmm.consensus) # TODO write to file and emit as output

    with open(tmp_hmm_path, "wb") as output_file:
        hmm.write(output_file)

    log_time(start_time, "run_hmmbuild (pyhmmer): ")

def mask_sequence(sequence, env_from, env_to):
    # Unpack the tuple (header, seq)
    header, seq = sequence[0]  # sequence[0] is a tuple

    # Modify header and sequence
    header = f"{header}/{env_from}_{env_to}"
    seq = seq[env_from - 1:env_to]  # -1 for 0-based indexing

    # Replace the old tuple with the new one
    sequence[0] = (header, seq)

    return sequence

def run_hmmsearch(pyhmmer_seqs, pyfastx_obj, exit_flag):
    """Runs HMMER's hmmsearch using pyhmmer."""
    start_time = time.time()

    filtered_sequences = []
    os.remove(tmp_family_sequences_path) # emptying initial MSA tmp file

    with pyhmmer.plan7.HMMFile(tmp_hmm_path) as hmm_file:
        hmm = hmm_file.read()
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

def run_hmmalign():
    """Runs HMMER's hmmalign using pyhmmer."""
    start_time = time.time()

    num_seqs_result = 0

    with pyhmmer.plan7.HMMFile(tmp_hmm_path) as hmm_file:
        hmm = hmm_file.read()
        with pyhmmer.easel.SequenceFile(tmp_family_sequences_path, digital=True) as seq_file:
            seqs = seq_file.read_block()
            hmmalign_res = pyhmmer.hmmer.hmmalign(hmm, seqs, trim=True)
            with open(tmp_align_msa_path, "wb") as outfile:
                hmmalign_res.write(outfile, format="stockholm") # TODO a2m here? (expected 'stockholm', 'pfam', 'a2m', 'psiblast', 'selex', 'afa', 'clustal', 'clustallike', 'phylip' or 'phylips')
                # TODO manipulate the object below
                # for name, aligned in zip(hmmalign_res.names, hmmalign_res.alignment):
                #     print(name.decode(), " ", aligned)
                num_seqs_result = len(hmmalign_res.names)

    log_time(start_time, "run_hmmalign (pyhmmer): ")

    return [name.decode() for name in hmmalign_res.names], num_seqs_result

def unmask_sequence_names(names):
    return [name.split('/')[0] for name in names]

def get_sequences_from_stockholm(file): # TODO remove and do with pyhmmer obj instead?
    with open(file, "r") as f:
        return AlignIO.read(f, "stockholm")

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

    shutil.copy(tmp_esl_weight_path, tmp_seed_msa_path) # TODO update with extra logic after Alex strategy's implemented

    log_time(start_time, "run_esl_weight: ")

def extract_RF():
    with open(tmp_seed_msa_path, 'r') as file: # TODO maybe change if .sto file is parsed or changed
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
    shutil.move(tmp_seed_msa_path, os.path.join(seed_msa_folder,  f'{arg_chunk_num}_{iteration}.sto')) # TODO tmp_seed_msa_path is currently tmp_seed_msa_sto_path
    shutil.move(tmp_align_msa_path,    os.path.join(align_msa_folder, f'{arg_chunk_num}_{iteration}.sto'))
    shutil.move(tmp_hmm_path,          os.path.join(hmm_folder,       f'{arg_chunk_num}_{iteration}.hmm'))
    shutil.move(tmp_domtblout_path,    os.path.join(domtblout_folder, f'{arg_chunk_num}_{iteration}.domtblout'))
    shutil.move(tmp_rf_path,           os.path.join(rf_folder,        f'{arg_chunk_num}_{iteration}.txt'))

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
        # TODO run trim_seed_msa / clipkit? Probably not needed
        # TODO renumber regions after clipping

        total_checked_sequences = family_members
        filtered_seq_names = []
        hand_flag = False
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
                run_hmmbuild(tmp_seed_msa_path)
                filtered_seq_names = run_hmmsearch(pyhmmer_seqs, mgnifams_pyfastx_obj, exit_flag)
                if (len(filtered_seq_names) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    break

                recruited_sequence_names, num_seqs_for_esl = run_hmmalign()
                unmasked_recruited_sequence_names = unmask_sequence_names(recruited_sequence_names)
                if (num_seqs_for_esl > 70000): # TODO remove if using mmseqs instead
                    discard_flag = True
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} too many sequences for esl ({num_seqs_for_esl}).\n")
                    break

                new_recruited_sequences = set(unmasked_recruited_sequence_names) - set(total_checked_sequences)
                if not new_recruited_sequences:
                    exit_flag = True
                    with open(log_file, 'a') as file:
                        file.write("Exiting-CONVERGED: no new sequences recruited.\n")
                    with open(converged_families_file, 'a') as file:
                        file.write(f"{iteration}\n")

            if exit_flag: # exit strategy branch
                with open(log_file, 'a') as file:
                    file.write("Exiting branch strategy:\n")
                run_hmmbuild(tmp_seed_msa_path, hand=hand_flag)
                extract_RF()
                filtered_seq_names = run_hmmsearch(pyhmmer_seqs, mgnifams_pyfastx_obj, exit_flag)
                if (len(filtered_seq_names) == 0): # low complexity sequence, confounding cluster, discard and move on to the next
                    discard_flag = True
                    break

                membership_percentage = check_seed_membership(original_sequence_names, unmask_sequence_names(filtered_seq_names))
                if (membership_percentage < 0.9): 
                    discard_flag = True
                    with open(log_file, 'a') as file:
                        file.write(f"Discard-Warning: {iteration} seed percentage in MSA is {membership_percentage}\n")
                    break
                elif (membership_percentage < 1):
                    with open(log_file, 'a') as file:
                        file.write(f"Warning: {iteration} seed percentage in MSA is {membership_percentage}\n")
                
                run_hmmalign() # final full MSA, including smaller sequences
                break

            # main strategy branch continue
            total_checked_sequences += list(new_recruited_sequences)
            with open(log_file, 'a') as file:
                file.write("total_checked_sequences calculated and starting run_esl_weight\n")
            run_esl_weight() # removes redundant sequences
            hand_flag = True # means the algorithm reached here at elast once, generating a stockholm format alignment

        # Exiting family loop
        if (discard_flag): # unsuccessfully
            with open(log_file, 'a') as file:
                file.write("Discarding cluster " + family_rep + "\n")
            
            with open(discarded_clusters_file, 'a') as outfile:  # TODO + reason for multiqc
                outfile.write(str(family_rep) + "\n")
            iteration -= 1 # keep proper track of family ids
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
