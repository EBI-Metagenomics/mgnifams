#!/usr/bin/env python

import sys
import gzip
from Bio import SeqIO

def get_matching_cluster_items(cluster_tsv_file, cluster_rep_id):
    with open(cluster_tsv_file) as f:
        cluster_items = []
        for line in f:
            columns = line.strip().split("\t")
            if columns[0] == cluster_rep_id:
                cluster_items.append(columns[1])
    return cluster_items

def filter_sequences(sequence_db_file, cluster_items):
    with open(sequence_db_file, "rt") as handle:
        filtered_records = []
        for record in SeqIO.parse(handle, "fasta"):
            is_uniprot = False
            # Check if the header is in Uniprot format
            if record.id.startswith('sp|'):
                seq_id = record.id.split('|')[1]  # Split on '|' and take the second element
                is_uniprot = True
            else:
                seq_id = record.id
            if seq_id in cluster_items:
                filtered_records.append(record)
                if is_uniprot:
                    global contains_sp
                    contains_sp = True
    return filtered_records

# Get command line arguments
cluster_tsv_file = sys.argv[1]
sequence_db_file = sys.argv[2]
output_file = sys.argv[3]
cluster_rep_id = sys.argv[4]
contains_sp = False

cluster_items = get_matching_cluster_items(cluster_tsv_file, cluster_rep_id)
filtered_records = filter_sequences(sequence_db_file, cluster_items)

if not contains_sp:
    output_file = "unknown/" + output_file
else:
    output_file = "known/" + output_file

with open(output_file, "w") as out_handle:
    SeqIO.write(filtered_records, out_handle, "fasta")