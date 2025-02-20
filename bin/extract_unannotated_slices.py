#!/usr/bin/env python3

import json
import csv
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Extract unannotated protein slices from input CSV.")
    parser.add_argument(
        "-i", "--input_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Path to the input CSV file containing protein sequences and metadata."
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Path to the output FASTA file to store the extracted sequences."
    )
    parser.add_argument(
        "-l", "--min_sequence_length",
        required=True,
        metavar="INT",
        type=int,
        help="Minimum sequence length for unannotated slices."
    )
    return parser.parse_args(args)

def hasAnnotation(metadata):
    return "p" in metadata or "m" in metadata

def sliceProtein(row, min_sequence_length):
    metadata, sequence, mgyp = row["metadata"], row["sequence"], row["mgyp"]
    pfams = metadata.get("p", [])
    mgnifams = metadata.get("m", [])

    # Create a sorted list of non-overlapping regions from both "p" and "m"
    sorted_regions = sorted((region[-2], region[-1]) for region in pfams + mgnifams)
    merged_regions = []
    
    for start, end in sorted_regions:
        if not merged_regions or start > merged_regions[-1][1]:
            merged_regions.append([start, end])
        else:
            merged_regions[-1][1] = max(merged_regions[-1][1], end)

    sliced_sequences = {}
    current_start = 1
    for start, end in merged_regions:
        if start - current_start >= min_sequence_length:
            key = f"{mgyp}_{current_start}_{start - 1}"
            sliced_sequences[key] = sequence[current_start - 1:start - 1]
        current_start = end + 1

    if len(sequence) - current_start + 1 >= min_sequence_length:
        key = f"{mgyp}_{current_start}_{len(sequence)}"
        sliced_sequences[key] = sequence[current_start - 1:]

    return sliced_sequences

def writeFastaSingleLine(records, file_handle):
    for record in records:
        file_handle.write(f">{record.id}\n{str(record.seq)}\n")

def main():
    args = parse_args()
    csv.field_size_limit(500000)
    
    with open(args.input_file, "r") as infile, open(args.output_file, "w") as outfile:
        csv_reader = csv.DictReader(infile)
        records_batch = []

        for row in csv_reader:
            mgyp = row["mgyp"]
            sequence = row["sequence"]
            metadata = row["metadata"]

            try:
                metadata = json.loads(metadata)  # Parse JSON once here
            except json.JSONDecodeError:
                print(f"Skipping row due to JSON error: {row}")
                continue
            
            if hasAnnotation(metadata):  # Pass parsed dict instead of JSON string
                row_dict = {"mgyp": mgyp, "sequence": sequence, "metadata": metadata}
                sliced_sequences = sliceProtein(row_dict, args.min_sequence_length)
                for key, seq in sliced_sequences.items():
                    records_batch.append(SeqRecord(Seq(seq), id=key, description=""))
            else:
                if len(sequence) >= args.min_sequence_length:
                    records_batch.append(SeqRecord(Seq(sequence), id=mgyp, description=""))

            if len(records_batch) >= 10000:  # Adjust batch size as needed
                writeFastaSingleLine(records_batch, outfile)
                records_batch.clear()

        # Write any remaining records
        if records_batch:
            writeFastaSingleLine(records_batch, outfile)

if __name__ == "__main__":
    main()
