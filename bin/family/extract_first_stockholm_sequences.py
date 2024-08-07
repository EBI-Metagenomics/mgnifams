import os
import sys
from Bio import AlignIO
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_file(file_path, output_fasta, problematic_seqs_file):
    basename = os.path.basename(file_path).split('.')[0]

    try:
        alignment = AlignIO.read(file_path, "stockholm")
        first_record = alignment[0]
        sequence = ''.join([c.upper() for c in first_record.seq if c.isalpha()])

        with open(output_fasta, 'a') as out_f:
            out_f.write(f">{basename}\n{sequence}\n")
        logging.info(f"Processed {file_path} successfully.")
    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}")
        with open(problematic_seqs_file, 'a') as out_f:
            out_f.write(f">{basename}\n")

def main(input_msa_folder, output_fasta, problematic_seqs_file):
    for file in os.listdir(input_msa_folder):
        process_file(os.path.join(input_msa_folder, file), output_fasta, problematic_seqs_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_first_stockholm_sequences.py <input_msa_folder> <output_fasta> <problematic_ids.txt>")
        sys.exit(1)

    with open(sys.argv[3], 'a'): # initialise non-optional problematic_ids.txt
        os.utime(sys.argv[3], None)
    main(sys.argv[1], sys.argv[2], sys.argv[3])