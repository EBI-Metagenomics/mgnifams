#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

def extract_first(file_path):
    with open(file_path) as f:
        for record in SeqIO.parse(f, "fasta"):
            return(record.id)

def process_file(file_path, output_mapping):
    basename            = os.path.basename(file_path).split('.')[0]
    first_sequence_name = extract_first(file_path)
    
    with open(output_mapping, 'a') as out_f:
        out_f.write(f"{basename},{first_sequence_name}\n")

def main(input_msa_folder, output_mapping):
    for file in os.listdir(input_msa_folder):
        process_file(os.path.join(input_msa_folder, file), output_mapping)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_first_stockholm_sequences.py <input_msa_folder> <output_mapping>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
