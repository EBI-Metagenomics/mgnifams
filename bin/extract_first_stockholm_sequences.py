import os
import sys
from Bio import AlignIO

def read_ids(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]

def process_file(file_path, ids, output_fasta, mode):
    basename = os.path.basename(file_path).split('.')[0]
    if mode == "unannotated" and basename not in ids:
        return

    alignment = AlignIO.read(file_path, "stockholm")
    first_record = alignment[0]
    sequence = ''.join([c.upper() for c in first_record.seq if c.isalpha()])

    with open(output_fasta, 'a') as out_f:
        out_f.write(f">{basename}_{first_record.id}\n{sequence}\n")

def main(input_msa_folder, unannotated_ids_path, mode, output_fasta):
    if mode not in ["unannotated", "all"]:
        raise ValueError("Invalid mode. Use 'unannotated' or 'all'.")
    
    ids = read_ids(unannotated_ids_path) if mode == "unannotated" else []
    for file in os.listdir(input_msa_folder):
        process_file(os.path.join(input_msa_folder, file), ids, output_fasta, mode)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_first_stockholm_sequences.py <input_msa_folder> <unannotated_ids_path> <mode> <output_fasta>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
