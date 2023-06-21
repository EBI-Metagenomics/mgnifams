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

# Load the gzipped sequence database
with gzip.open(sequence_db_file, "rt") as handle:
    seq_db = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

# Write the cluster sequences to a new FASTA file
with open(output_file, "w") as f:
    for seq_id in cluster_seqs:
        SeqIO.write(seq_db[seq_id], f, "fasta")
