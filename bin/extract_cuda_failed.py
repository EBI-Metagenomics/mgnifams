#!/usr/bin/env python3

import os
import argparse

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

def extract_fasta_sequences(input_fasta, failed_sequences, output_file):
    with open(input_fasta, 'r') as infile, open(output_file, 'w') as outfile:
        write_sequence = False
        for line in infile:
            if line.startswith('>'):
                protein_id = line.split()[0][1:].replace("/", "_")  # Extract the protein ID and replace / with _
                if protein_id in failed_sequences:
                    write_sequence = True
                    outfile.write(line)  # Write the header line
                else:
                    write_sequence = False
            elif write_sequence:
                outfile.write(line)  # Write the sequence lines

def main():
    parser = argparse.ArgumentParser(description="Extract sequences that failed due to CUDA memory errors.")
    parser.add_argument("--input_fasta", help="Input FASTA file")
    parser.add_argument("--input_scores_folder", help="Folder containing logs with CUDA memory errors")
    parser.add_argument("--output_fasta", help="Output FASTA file with failed sequences")
    args = parser.parse_args()

    failed_sequences = extract_failed_sequences(args.input_scores_folder)
    extract_fasta_sequences(args.input_fasta, failed_sequences, args.output_fasta)

if __name__ == "__main__":
    main()
