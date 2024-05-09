import os
import sys
from Bio import AlignIO

def process_file(file_path, output_fasta):
    basename = os.path.basename(file_path).split('.')[0]

    alignment = AlignIO.read(file_path, "stockholm")
    first_record = alignment[0]
    sequence = ''.join([c.upper() for c in first_record.seq if c.isalpha()])

    with open(output_fasta, 'a') as out_f:
        out_f.write(f">{basename}\n{sequence}\n")

def main(input_msa_folder, output_fasta):
    for file in os.listdir(input_msa_folder):
        process_file(os.path.join(input_msa_folder, file), output_fasta)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_first_stockholm_sequences.py <input_msa_folder> <output_fasta>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
