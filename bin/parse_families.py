import os
import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    print("Usage: python3 parse_families.py <fastaFile> <clust_tsv>")
    sys.exit(1)

fasta_file_path = sys.argv[1]
cluster_file_path = sys.argv[2]

# Load the fasta file into a dictionary with SeqIO (Biopython)
fasta_records = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file_path, "fasta")}

# Initialize
current_rep = None
current_cluster = []
prev_rep = None  # Variable to keep track of the previously encountered representative

# Ensure the non_singletons directory exists
output_dir = "non_singletons"
os.makedirs(output_dir, exist_ok=True)

# Open files for representatives and prepare for writing
with open("reps.fa", 'w') as reps_f, open("rep_names.txt", 'w') as rep_names_f, open("singletons.txt", 'w') as singletons_f:
    # Read cluster file line by line
    with open(cluster_file_path, 'r') as f:
        for line in f:
            rep, member = line.strip().split('\t')
            
            # Check if we moved to a new cluster
            if current_rep and rep != current_rep:
                if len(current_cluster) == 1:  # Singleton cluster
                    singletons_f.write(current_rep + "\n")
                else:  # Non-singleton cluster
                    with open(os.path.join(output_dir, current_rep + ".fasta"), 'w') as f_out:
                        for member_id in current_cluster:
                            f_out.write(">" + member_id + "\n")
                            f_out.write(fasta_records[member_id] + "\n")
                current_cluster = []

            # Update variables
            if rep != prev_rep:  # Check if rep has not been written yet
                reps_f.write(">" + rep + "\n")
                reps_f.write(fasta_records[rep] + "\n")
                rep_names_f.write(rep + "\n")
                prev_rep = rep  # Update prev_rep to avoid writing duplicates
            
            current_rep = rep
            current_cluster.append(member)
            
        # Handle the last cluster
        if len(current_cluster) == 1:
            singletons_f.write(current_rep + "\n")
        else:
            with open(os.path.join(output_dir, current_rep + ".fasta"), 'w') as f_out:
                for member_id in current_cluster:
                    f_out.write(">" + member_id + "\n")
                    f_out.write(fasta_records[member_id] + "\n")
