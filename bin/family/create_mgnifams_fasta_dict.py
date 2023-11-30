import sys
import dill
from Bio import SeqIO

import time # benchmarking, TODO remove

def parse_args():
    if not (len(sys.argv) == 2):
        print("Usage: python3 create_mgnifams_fasta_dict.py [MGnifams FASTA file]")
        sys.exit(1)

    globals().update({
        "mgnifams_input_file" : sys.argv[1]
    })      

def define_globals():
    globals().update({
        "log_file"               : "log.txt",
        "mgnifams_dict_pkl_file" : "mgnifams_dict.pkl"
    })
    
def create_mgnifams_fasta_dict():
    start_time = time.time()

    mgnifams_fasta_dict = {record.id: record for record in SeqIO.parse(mgnifams_input_file, "fasta")}
    with open(mgnifams_dict_pkl_file, 'wb') as file:
        dill.dump(mgnifams_fasta_dict, file)

    with open(log_file, 'a') as file:
        file.write("create_mgnifams_fasta_dict: ")
        file.write(str(time.time() - start_time) + "\n")

    return mgnifams_fasta_dict

def main():
    parse_args()
    define_globals()

    mgnifams_fasta_dict = create_mgnifams_fasta_dict()

if __name__ == "__main__":
    main()
