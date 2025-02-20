#!/usr/bin/env python3

import json
import csv
import argparse

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Extract unannotated protein slices from input CSV.")
    parser.add_argument("-i", "--input_file", required=True, type=str, help="Input CSV file.")
    parser.add_argument("-o", "--output_file", required=True, type=str, help="Output FASTA file.")
    parser.add_argument("-l", "--min_sequence_length", required=True, type=int, help="Min sequence length.")
    return parser.parse_args(args)

def sliceProtein(mgyp, sequence, metadata, min_sequence_length):
    pfams = metadata.get("p", [])
    mgnifams = metadata.get("m", [])

    sorted_regions = sorted((region[-2], region[-1]) for region in pfams + mgnifams)
    merged_regions = []
    
    for start, end in sorted_regions:
        if not merged_regions or start > merged_regions[-1][1]:
            merged_regions.append([start, end])
        else:
            merged_regions[-1][1] = max(merged_regions[-1][1], end)

    sliced_sequences = []
    current_start = 1
    for start, end in merged_regions:
        if start - current_start >= min_sequence_length:
            sliced_sequences.append(f">{mgyp}_{current_start}_{start - 1}\n{sequence[current_start - 1:start - 1]}\n")
        current_start = end + 1

    if len(sequence) - current_start + 1 >= min_sequence_length:
        sliced_sequences.append(f">{mgyp}_{current_start}_{len(sequence)}\n{sequence[current_start - 1:]}\n")

    return sliced_sequences

def main():
    args = parse_args()
    csv.field_size_limit(500000)
    
    with open(args.input_file, "r") as infile, open(args.output_file, "w") as outfile:
        csv_reader = csv.DictReader(infile)
        buffer = []

        for row in csv_reader:
            try:
                metadata = json.loads(row["metadata"])
            except json.JSONDecodeError:
                print(f"Skipping row due to JSON error: {row}")
                continue

            sequence, mgyp = row["sequence"], row["mgyp"]
            if "p" in metadata or "m" in metadata:
                buffer.extend(sliceProtein(mgyp, sequence, metadata, args.min_sequence_length))
            elif len(sequence) >= args.min_sequence_length:
                buffer.append(f">{mgyp}\n{sequence}\n")

            if len(buffer) >= 10000:
                outfile.writelines(buffer)
                buffer.clear()

        if buffer:
            outfile.writelines(buffer)

if __name__ == "__main__":
    main()
