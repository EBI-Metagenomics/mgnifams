#!/usr/bin/env python

import sys
from Bio import SeqIO

# Get command line arguments
fasta_file = sys.argv[1]
output_file = sys.argv[2]
unknown_flag = bool(int(sys.argv[3]))
contains_sp = False

# Load the FASTA file
with open(fasta_file, "rt") as handle:
    records = list(SeqIO.parse(handle, "fasta"))
    for record in records:
        # Check if the sequence ID starts with 'sp'
        if record.id.startswith('sp'):
            contains_sp = True
            break

if unknown_flag:
    # If no sequence starts with 'sp', write the file contents to output_file
    if not contains_sp:
        with open(output_file, "w") as out_handle:
            for record in records:
                SeqIO.write(record, out_handle, "fasta")
else:
    if contains_sp:
        with open(output_file, "w") as out_handle:
            for record in records:
                SeqIO.write(record, out_handle, "fasta")
