#!/usr/bin/env python

import sys
import gzip
from Bio import SeqIO

# Get command line arguments
cluster_tsv_file = sys.argv[1]
sequence_db_file = sys.argv[2]
output_file = sys.argv[3]
cluster_rep_id = sys.argv[4]

# Load the cluster TSV file
with open(cluster_tsv_file) as f:
    clusters = {}
    for line in f:
        rep_id, seq_id = line.strip().split("\t")
        if rep_id not in clusters:
            clusters[rep_id] = []
        clusters[rep_id].append(seq_id)

# Extract the sequences for a specific cluster
cluster_seqs = clusters[cluster_rep_id]

# Load the unzipped sequence database
with open(sequence_db_file, "rt") as handle: # if zipped: gzip.open instead
    seq_db = {}
    for record in SeqIO.parse(handle, "fasta"):
        # Check if the header is in Uniprot format
        if record.id.startswith('sp|'):
            seq_id = record.id.split('|')[1]  # Split on '|' and take the second element
        else:
            seq_id = record.id
        seq_db[seq_id] = record

# Write the cluster sequences to a new FASTA file
with open(output_file, "w") as f:
    for seq_id in cluster_seqs:
        SeqIO.write(seq_db[seq_id], f, "fasta")
