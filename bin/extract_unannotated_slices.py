#!/usr/bin/env python3

import json
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def hasAnnotation(metadata):
    try:
        data = json.loads(metadata)
        return "p" in data or "m" in data
    except json.JSONDecodeError:
        return False

def sliceProtein(row, min_sequence_length):
    metadata, sequence, mgyp = row["metadata"], row["sequence"], row["mgyp"]
    data = json.loads(metadata)
    
    # Extract Pfam ("p") and "m" regions
    pfams = data.get("p", [])
    mgnifams = data.get("m", [])
    
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

def main(input_file, output_file, min_sequence_length):
    csv.field_size_limit(500000)
    
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        csv_reader = csv.reader(infile)
        next(csv_reader, None) # Skip the header
        records_batch = []

        for row in csv_reader:
            if len(row) != 5:
                print(f"Skipping malformed line: {row}")
                continue
            mgyp, sequence, _, _, metadata = row
            if hasAnnotation(metadata):
                row_dict = {"mgyp": mgyp, "sequence": sequence, "metadata": metadata}
                sliced_sequences = sliceProtein(row_dict, min_sequence_length)
                for key, seq in sliced_sequences.items():
                    records_batch.append(SeqRecord(Seq(seq), id=key, description=""))
            else:
                if len(sequence) >= min_sequence_length:
                    records_batch.append(SeqRecord(Seq(sequence), id=mgyp, description=""))

            if len(records_batch) >= 10000:  # Adjust batch size as needed
                writeFastaSingleLine(records_batch, outfile)
                records_batch.clear()

        # Write any remaining records
        if records_batch:
            writeFastaSingleLine(records_batch, outfile)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_unannotated_slices.py <input_file> <output_file> <min_sequence_length>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
