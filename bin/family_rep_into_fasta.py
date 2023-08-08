#!/usr/bin/env python

import sys
from Bio import SeqIO

def get_matching_cluster_ids(clusters_tsv, cluster_rep_id):
    with open(clusters_tsv) as f:
        cluster_ids = set()
        found_cluster = False
        for line in f:
            columns = line.strip().split("\t")
            if columns[0] == cluster_rep_id:
                cluster_ids.add(columns[1])
                found_cluster = True
            # Break if we've moved past the target cluster
            elif found_cluster:
                break
    
    return cluster_ids

def find_cluster_sequences(input_fasta, cluster_ids):
    contains_sp = False
    cluster_seqs = []

    with open(input_fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id.split('|')[1] if record.id.startswith('sp|') else record.id
            if seq_id in cluster_ids:
                cluster_seqs.append(record)
                if record.id.startswith('sp|'):
                    contains_sp = True
                cluster_ids.remove(seq_id) # Remove the found item from the remaining set
                if not cluster_ids: # If no more items are remaining, break
                    break
    
    return contains_sp, cluster_seqs

def write_cluster_sequences(contains_sp, cluster_seqs, output_file):
    output_file = (contains_sp and "known/" or "unknown/") + output_file

    with open(output_file, "w") as out_handle:
        SeqIO.write(cluster_seqs, out_handle, "fasta")


# Get command line arguments
clusters_tsv = sys.argv[1]
input_fasta = sys.argv[2]
output_file = sys.argv[3]
cluster_rep_id = sys.argv[4]

cluster_ids = get_matching_cluster_ids(clusters_tsv, cluster_rep_id)
contains_sp, cluster_seqs = find_cluster_sequences(input_fasta, cluster_ids)
write_cluster_sequences(contains_sp, cluster_seqs, output_file)
