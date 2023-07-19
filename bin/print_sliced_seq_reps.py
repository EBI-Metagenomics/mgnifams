#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def read_slices(slices_file):
    col_names = ['from', 'to']
    slices = pd.read_csv(slices_file, header=None, names=col_names)
    return(slices)

def get_first_sequence(msa_file):
    with open(msa_file, 'r') as file:
        sequences = file.read().split('>')[1:]  # skip the empty string before the first '>'
        first_sequence = sequences[0].split('\n', 1)[1]  # get the sequence part of the first sequence
        first_sequence = first_sequence.replace('\n', '')  # remove newlines
        return first_sequence

def get_substrings(sequence, slices):
    substrings = [sequence[row['from']-1 : row['to']] for index, row in slices.iterrows()]
    return substrings

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scan a sliced mafft result first sequence for domains.')
    parser.add_argument('slices_file', type=str, help='The path to the slice ranges')
    parser.add_argument('msa_file', type=str, help='The path to the family msa')
    parser.add_argument('output_file', type=str, help='The path to the output file')
    parser.add_argument('min_slice_length', type=int, help='The minimum length of unannotated slices to keep')

    args = parser.parse_args()

    slices = read_slices(args.slices_file)
    sequence = get_first_sequence(args.msa_file)
    substrings = get_substrings(sequence, slices)
    msa_filename = os.path.basename(args.msa_file)

    
    with open(args.output_file, 'w') as f:
        # Iterate over substrings and their corresponding slices
        for s, (index, row) in zip(substrings, slices.iterrows()):
            if len(s) >= args.min_slice_length:
                f.write(f">{msa_filename}_{row['from']}_{row['to']}\n{s}\n")

    if os.path.getsize(args.output_file) == 0: # remove file if nothing written, optional NF output
        os.remove(args.output_file)