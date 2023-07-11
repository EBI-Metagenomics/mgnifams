#!/usr/bin/env python3

import pandas as pd
import argparse

def load_file(file_path):
    # specify column names
    column_names = ['target_name', 'target_accession', 'tlen', 'query_name', 'query_accession', 'qlen', 'E-value', 
                    'seq_score', 'seq_bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'dom_score', 'dom_bias', 'hmm_from', 'hmm_to', 
                    'ali_from', 'ali_to', 'env_from', 'env_to', 'acc', 'description_of_target']

    # read the file
    df = pd.read_csv(file_path, delim_whitespace=True, comment='#', names=column_names)
    df = df.sort_values(by=['hmm_from'])

    # print the dataframe
    return(df)

def init_merged_list(df):
    # List to hold the merged intervals
    merged = []

    for row in df.itertuples():
        if not merged or merged[-1][1] < row.hmm_from:
            # if the list of merged intervals is empty, or the current interval does not overlap with the previous, append it
            merged.append([row.hmm_from, row.hmm_to])
        else:
            # otherwise, there is overlap, so we merge the current and previous intervals
            merged[-1][1] = max(merged[-1][1], row.hmm_to)
    
    return(merged)

# Create a function to get the uncovered ranges
def find_uncovered(covered, total_length):
    uncovered = []
    last_end = 0

    for start, end in covered:
        if start > last_end + 1:
            uncovered.append((last_end + 1, start - 1))
        last_end = end

    if last_end < total_length:
        uncovered.append((last_end + 1, total_length))

    return uncovered

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Load a HMMER output file into a pandas DataFrame.')
    parser.add_argument('file_path', type=str, help='The path to the HMMER output file')
    parser.add_argument('output_file', type=str, help='The path to the output file')

    args = parser.parse_args()
    df = load_file(args.file_path)
    merged = init_merged_list(df)

    total_length = df['tlen'].iloc[0]
    uncovered = find_uncovered(merged, total_length)

    with open(args.output_file, 'w') as f:
        for interval in uncovered:
            f.write(f'{interval[0]},{interval[1]}\n')