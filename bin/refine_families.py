import sys
import os
import subprocess
import pandas as pd
from Bio import SeqIO

import time # benchmarking, TODO remove

def create_family_dataframe(families_tsv):
    df = pd.read_csv(families_tsv, sep='\t', header=None, names=['representative', 'member'])
    family_sizes = df.groupby('representative').size()
    bookkeeping_df = df.set_index('representative').join(family_sizes.rename('size'), on='representative')

    return bookkeeping_df

def find_next_largest_family(bookkeeping_df):
    unchecked_families = bookkeeping_df[~bookkeeping_df['checked']]
    if unchecked_families.empty:
        return None, None

    largest_family_rep = unchecked_families['size'].idxmax()
    largest_family_members = unchecked_families.loc[largest_family_rep]['member']
    bookkeeping_df.at[largest_family_rep, 'checked'] = True

    return largest_family_rep, largest_family_members.tolist() if isinstance(largest_family_members, pd.Series) else [largest_family_members]

def get_family_fasta(members, fasta_dict, output_fasta):
    sequences_to_write = [fasta_dict[member] for member in members if member in fasta_dict]

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(sequences_to_write, output_handle, "fasta")

def generate_msa(input_fasta, output_msa):
    mafft_command = ["mafft", "--quiet", "--auto", input_fasta]
    with open(output_msa, "w") as output_handle:
        subprocess.run(mafft_command, stdout=output_handle)

def generate_hmm(msa_file, output_hmm):
    hmmbuild_command = ["hmmbuild", "--informat", "afa", output_hmm, msa_file]
    subprocess.run(hmmbuild_command, stdout=subprocess.DEVNULL)

def recruit_sequences(hmm_file, fasta_file, output_recruitment):
    hmmsearch_command = ["hmmsearch", "--tblout", output_recruitment, hmm_file, fasta_file]
    subprocess.run(hmmsearch_command, stdout=subprocess.DEVNULL)

def filter_recruited(recruitment_file, evalue_threshold, bitscore_threshold, checked_sequences):
    filtered_sequences = []

    with open(recruitment_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split()
                evalue = float(columns[4])
                bitscore = float(columns[5])
                sequence_name = columns[0]

                if sequence_name not in checked_sequences and evalue < evalue_threshold and bitscore >= bitscore_threshold:
                    filtered_sequences.append(sequence_name)

    return filtered_sequences

# TODO try for larger inputs
# def filter_recruited(recruitment_file, evalue_threshold, bitscore_threshold, checked_sequences):
#     df = pd.read_csv(recruitment_file, sep='\s+', comment='#', header=None, usecols=[0, 4, 5])
#     df.columns = ['sequence_name', 'evalue', 'bitscore']

#     filtered_df = df[(~df['sequence_name'].isin(checked_sequences)) & 
#                      (df['evalue'] < evalue_threshold) & 
#                      (df['bitscore'] >= bitscore_threshold)]

#     return filtered_df['sequence_name'].tolist()

def write_to_file(output_file, family_rep, family_members):
    lines = [f"{family_rep}\t{member}\n" for member in family_members]
    with open(output_file, 'a') as file:
        file.writelines(lines)

def update_bookkeeping(bookkeeping_df, largest_family):
    # Vectorized update for 'checked' status
    bookkeeping_df.loc[bookkeeping_df['member'].isin(largest_family), 'checked'] = True
    # Efficiently recalculating family sizes
    bookkeeping_df['size'] = bookkeeping_df.groupby(level=0)['checked'].apply(lambda x: (~x).sum())

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 refine_families.py [Families TSV] [FASTA file] [Minimum number of family members] [Output file]")
        sys.exit(1)

    families_tsv = sys.argv[1]
    fasta_file = sys.argv[2]
    minimum_members = int(sys.argv[3])
    output_file = sys.argv[4]

    # Benchmarking # TODO remove
    # total_time_reading_input_fasta = 0
    # total_time_get_family_fasta = 0
    # total_time_generate_msa = 0
    # total_time_generate_hmm = 0
    # total_time_recruit_sequences = 0
    # total_time_filter_recruited = 0
    # total_time_write_to_file = 0
    # total_time_update_bookkeeping = 0
    ###

    tmp_folder = "tmp"
    msa_folder = "msa"
    hmm_folder = "hmm"
    for folder in [tmp_folder, msa_folder, hmm_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
        if not os.path.exists(tmp_folder):
            os.makedirs(tmp_folder)
    family_sequences_path = os.path.join(tmp_folder, 'family_sequences.fasta')
    family_alignment_path = os.path.join(tmp_folder, 'family_alignment.msa')
    family_model_path = os.path.join(tmp_folder, 'family_model.hmm')
    recruited_sequences_path = os.path.join(tmp_folder, 'recruited_sequences.txt')

    evalue_threshold = 0.001
    bitscore_threshold = 20

    bookkeeping_df = create_family_dataframe(families_tsv)
    bookkeeping_df['checked'] = False
    checked_sequences = []

    # start_time = time.time()
    # Loading mgnifams_input.fa in memory, to reduce I/O in every iteration
    fasta_dict = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
    # total_time_reading_input_fasta += time.time() - start_time

    while True:
        family_rep, largest_family = find_next_largest_family(bookkeeping_df)
        if not largest_family or len(largest_family) < minimum_members:
            break

        print(f"Processing family: {family_rep}")
        inner_loop_counter = 0
        while True:
            inner_loop_counter += 1
            print(inner_loop_counter)
            
            # start_time = time.time()
            get_family_fasta(largest_family, fasta_dict, family_sequences_path)
            # total_time_get_family_fasta += time.time() - start_time

            # start_time = time.time()
            generate_msa(family_sequences_path, family_alignment_path)
            # total_time_generate_msa += time.time() - start_time

            # start_time = time.time()
            generate_hmm(family_alignment_path, family_model_path)
            # total_time_generate_hmm += time.time() - start_time

            # start_time = time.time()
            recruit_sequences(family_model_path, fasta_file, recruited_sequences_path)
            # total_time_recruit_sequences += time.time() - start_time

            # start_time = time.time()
            filtered_seq_names = filter_recruited(recruited_sequences_path, evalue_threshold, bitscore_threshold, checked_sequences)
            # total_time_filter_recruited += time.time() - start_time

            new_sequences = set(filtered_seq_names) - set(largest_family)
            if new_sequences:
                largest_family.extend(new_sequences)
                continue
            else:
                # start_time = time.time()
                write_to_file(output_file, family_rep, largest_family)
                # total_time_write_to_file += time.time() - start_time

                checked_sequences.extend(largest_family) # TODO keep this? not hcecking for already family recruited sequences in next loops

                # start_time = time.time()
                update_bookkeeping(bookkeeping_df, largest_family) # flag as checked and recalculate sizes in bookkeeping
                # total_time_update_bookkeeping += time.time() - start_time
                # TODO keep msa and hmm? better not at this point
                # os.rename(family_alignment_path, os.path.join(msa_folder, f"{family_rep}.msa"))
                # os.rename(family_model_path, os.path.join(hmm_folder, f"{family_rep}.hmm"))
                break

    # print(f"Total time for get_family_fasta: {# total_time_get_family_fasta} seconds")
    # print(f"Total time for generate_msa: {# total_time_generate_msa} seconds")
    # print(f"Total time for generate_hmm: {# total_time_generate_hmm} seconds")
    # print(f"Total time for recruit_sequences: {# total_time_recruit_sequences} seconds")
    # print(f"Total time for filter_recruited: {# total_time_filter_recruited} seconds")
    # print(f"Total time for write_to_file: {# total_time_write_to_file} seconds")
    # print(f"Total time for update_bookkeeping: {# total_time_update_bookkeeping} seconds")

if __name__ == "__main__":
    main()
