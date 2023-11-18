import bz2
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import csv

def hasPfam(metadata):
    try:
        data = json.loads(metadata)
        return "p" in data
    except json.JSONDecodeError:
        return False

def sliceProtein(row, min_slice_length):
    metadata, sequence, mgyp = row["metadata"], row["sequence"], row["mgyp"]
    data = json.loads(metadata)
    pfams = data.get("p", [])
    
    pfams.sort(key=lambda x: x[-2])

    sliced_sequences = {}
    current_start = 1
    for pfam in pfams:
        start, end = pfam[-2], pfam[-1]
        if start - current_start >= min_slice_length:
            key = f"{mgyp}_{current_start}_{start - 1}"
            sliced_sequences[key] = sequence[current_start - 1:start - 1]
        current_start = end + 1

    # Check the remaining sequence after the last pfam
    if len(sequence) - current_start + 1 >= min_slice_length:
        key = f"{mgyp}_{current_start}_{len(sequence)}"
        sliced_sequences[key] = sequence[current_start - 1:]

    return sliced_sequences

def writeFastaSingleLine(record, file):
    with open(file, 'a') as f:
        f.write(f">{record.id}\n{str(record.seq)}\n")

def main(input_file, output_file, min_slice_length):
    min_slice_length = int(min_slice_length)
    with bz2.open(input_file, "rt") as infile, open(output_file, "w") as outfile:
        csv_reader = csv.reader(infile)
        next(csv_reader, None) # Skip the header line
        for row in csv_reader:
            if len(row) != 5:
                print(f"Skipping malformed line: {row}")
                continue
            mgyp, sequence, _, _, metadata = row  # Unpack the row as a list
            if hasPfam(metadata):
                row_dict = {"mgyp": mgyp, "sequence": sequence, "metadata": metadata}
                sliced_sequences = sliceProtein(row_dict, min_slice_length)
                for key, seq in sliced_sequences.items():
                    record = SeqRecord(Seq(seq), id=key, description="")
                    writeFastaSingleLine(record, output_file)
            else:
                record = SeqRecord(Seq(sequence), id=mgyp, description="")
                writeFastaSingleLine(record, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 filter_unannotated_slices_fasta.py <input_file> <output_file> <min_slice_length>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
