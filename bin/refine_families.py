import sys
import os
import subprocess
import pandas as pd
from Bio import SeqIO

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

def get_family_fasta(members, fasta_file, output_fasta):
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in members:
                SeqIO.write(record, output_handle, "fasta")

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

def write_to_file(output_file, family_rep, family_members):
    with open(output_file, 'a') as file:
        for member in family_members:
            file.write(f"{family_rep}\t{member}\n")

def update_bookkeeping(bookkeeping_df, largest_family):
    # Marking sequences as checked
    for sequence in largest_family:
        if sequence in bookkeeping_df.index:
            bookkeeping_df.at[sequence, 'checked'] = True

    # Recalculating the family sizes
    bookkeeping_df['size'] = bookkeeping_df.groupby('representative')['checked'].transform(lambda x: (x == False).sum())

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 refine_families.py [Families TSV] [FASTA file] [Minimum number of family members] [Output file]")
        sys.exit(1)

    families_tsv = sys.argv[1]
    fasta_file = sys.argv[2]
    minimum_members = int(sys.argv[3])
    output_file = sys.argv[4]

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

    while True:
        family_rep, largest_family = find_next_largest_family(bookkeeping_df)
        if not largest_family or len(largest_family) < minimum_members:
            break

        print(f"Processing family: {family_rep}")
        inner_loop_counter = 0
        while True:
            inner_loop_counter += 1
            print(inner_loop_counter)
            
            get_family_fasta(largest_family, fasta_file, family_sequences_path)
            generate_msa(family_sequences_path, family_alignment_path)
            generate_hmm(family_alignment_path, family_model_path)
            recruit_sequences(family_model_path, fasta_file, recruited_sequences_path)
            filtered_seq_names = filter_recruited(recruited_sequences_path, evalue_threshold, bitscore_threshold, checked_sequences)

            new_sequences = set(filtered_seq_names) - set(largest_family)
            if new_sequences:
                largest_family.extend(new_sequences)
                continue
            else:
                write_to_file(output_file, family_rep, largest_family)
                checked_sequences.extend(largest_family) # TODO keep this? not hcecking for already family recruited sequences in next loops
                update_bookkeeping(bookkeeping_df, largest_family)# flag as checked and recalculate sizes in bookkeeping
                # TODO keep msa and hmm? better not at this point
                # os.rename(family_alignment_path, os.path.join(msa_folder, f"{family_rep}.msa"))
                # os.rename(family_model_path, os.path.join(hmm_folder, f"{family_rep}.hmm"))
                break

if __name__ == "__main__":
    main()
