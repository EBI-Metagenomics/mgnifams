import argparse
from Bio import SeqIO
import csv

def fasta_to_csv(fasta_file, csv_file):
    with open(csv_file, 'w', newline='') as csvfile:
        fieldnames = ['ID', 'Sequence']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = seq_record.id
            # Modify sequence ID if it starts with "sp|"
            if seq_id.startswith("sp|"):
                seq_id = seq_id.split("|")[1]

            writer.writerow({'ID': seq_id, 'Sequence': str(seq_record.seq)})

def main():
    parser = argparse.ArgumentParser(description='Convert a FASTA file to a CSV file.')
    parser.add_argument('fasta_file', help='The path to the input FASTA file.')
    parser.add_argument('csv_file', help='The path to the output CSV file.')
    args = parser.parse_args()

    fasta_to_csv(args.fasta_file, args.csv_file)

if __name__ == "__main__":
    main()
