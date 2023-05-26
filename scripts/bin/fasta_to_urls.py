#!/usr/bin/env python3

from Bio import SeqIO
import sys

fasta_file = sys.argv[1]
base_url = sys.argv[2]

with open(fasta_file) as fasta_file, open('urls.txt', 'w') as url_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id
        url = f"{base_url}/{protein_id}\n"
        url_file.write(url)
