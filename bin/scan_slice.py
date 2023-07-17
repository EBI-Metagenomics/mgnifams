#!/usr/bin/env python3

import argparse
import pandas as pd

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

    args = parser.parse_args()

    slices = read_slices(args.slices_file)
    sequence = get_first_sequence(args.msa_file)
    substrings = get_substrings(sequence, slices)
    # print(substrings)

    # TODO keep substrings over threshold characters, hmm and pyhmmer against uniprot