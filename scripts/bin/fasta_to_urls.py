#!/usr/bin/env python3

from Bio import SeqIO
import sys
import io

base_url = sys.argv[1]

fasta_content = sys.stdin.read()
fasta_file = io.StringIO(fasta_content)

with open('urls.txt', 'w') as url_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id
        url = f"{base_url}/{protein_id}\n"
        url_file.write(url)