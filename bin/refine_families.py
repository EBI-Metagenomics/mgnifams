import sys
import pandas as pd
from Bio import SeqIO
import subprocess

def create_family_dataframe(families_tsv):
    df = pd.read_csv(families_tsv, sep='\t', header=None, names=['representative', 'member'])
    family_sizes = df.groupby('representative').size()
    bookkeeping_df = df.set_index('representative').join(family_sizes.rename('size'), on='representative')

    return bookkeeping_df

def find_largest_family(bookkeeping_df):
    largest_family_rep = bookkeeping_df['size'].idxmax()
    largest_family_members = bookkeeping_df.loc[largest_family_rep]['member']

    return largest_family_members.tolist() if isinstance(largest_family_members, pd.Series) else [largest_family_members]

def get_family_fasta(members, fasta_file, output_fasta):
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in members:
                SeqIO.write(record, output_handle, "fasta")

def generate_msa(input_fasta, output_msa):
    mafft_command = ["mafft", "--auto", input_fasta]
    with open(output_msa, "w") as output_handle:
        subprocess.run(mafft_command, stdout=output_handle)

def generate_hmm(msa_file, output_hmm):
    hmmbuild_command = ["hmmbuild", "--informat", "afa", output_hmm, msa_file]
    subprocess.run(hmmbuild_command)

def recruit_sequences(hmm_file, fasta_file, output_recruitment):
    hmmsearch_command = ["hmmsearch", "--tblout", output_recruitment, hmm_file, fasta_file]
    subprocess.run(hmmsearch_command)

def filter_recruited(recruitment_file, evalue_threshold, bitscore_threshold):
    filtered_sequences = []

    with open(recruitment_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split()
                evalue = float(columns[4])
                bitscore = float(columns[5])

                if evalue <= evalue_threshold and bitscore >= bitscore_threshold:
                    sequence_name = columns[0]
                    filtered_sequences.append(sequence_name)

    return filtered_sequences

# TODO remove
def write_to_file(output_file, largest_family):
    with open(output_file, 'w') as file:
        file.write('\n'.join(largest_family))

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 refine_families.py [Families TSV] [FASTA file] [Minimum number of family members] [Output file]")
        sys.exit(1)

    families_tsv = sys.argv[1]
    fasta_file = sys.argv[2]
    minimum_members = int(sys.argv[3])
    output_file = sys.argv[4]

    evalue_threshold = 0.05
    bitscore_threshold = 20

    bookkeeping_df = create_family_dataframe(families_tsv)
    largest_family = find_largest_family(bookkeeping_df)
    write_to_file(output_file, largest_family) # TODO remove
    get_family_fasta(largest_family, fasta_file, 'family_sequences.fasta')
    generate_msa('family_sequences.fasta', 'family_alignment.msa')
    generate_hmm('family_alignment.msa', 'family_model.hmm')
    recruit_sequences('family_model.hmm', fasta_file, 'recruited_sequences.txt')
    filtered_seq_names = filter_recruited('recruited_sequences.txt', evalue_threshold, bitscore_threshold)
    get_family_fasta(filtered_seq_names, fasta_file, 'filtered_family_sequences.fasta')

if __name__ == "__main__":
    main()
