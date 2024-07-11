#!/usr/bin/env python3

import os
import sys

def extract_failed_sequences(input_folder):
    failed_sequences = []
    
    for filename in os.listdir(input_folder):
        file_path = os.path.join(input_folder, filename)
        if os.path.isfile(file_path):
            with open(file_path, 'r') as infile:
                for line in infile:
                    if "Failed (CUDA out of memory)" in line:
                        parts = line.split()
                        if "sequence" in parts:
                            sequence_index = parts.index("sequence") + 1
                            if sequence_index < len(parts):
                                failed_sequences.append(parts[sequence_index])
    
    return failed_sequences

def extract_fasta_sequences(input_folder, failed_sequences, output_file):
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_folder):
            file_path = os.path.join(input_folder, filename)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as infile:
                    write_sequence = False
                    for line in infile:
                        if line.startswith('>'):
                            protein_id = line.split()[0][1:]  # Extract the protein ID
                            if protein_id in failed_sequences:
                                write_sequence = True
                                outfile.write(line)  # Write the header line
                            else:
                                write_sequence = False
                        elif write_sequence:
                            outfile.write(line)  # Write the sequence lines

def main(input_fasta_folder, input_scores_folder, output_fasta):
    failed_sequences = extract_failed_sequences(input_scores_folder)
    extract_fasta_sequences(input_fasta_folder, failed_sequences, output_fasta)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_long_fas.py <input_fasta_folder> <input_scores_folder> <output_fasta>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
