import os
import sys
from Bio import AlignIO

def parse_protein_region(protein_id):
    number_of_underscores = protein_id.count('_')
    if (number_of_underscores == 0):
        protein = protein_id
        region  = "-"
    elif (number_of_underscores == 1):
        parts   = protein_id.split('/')
        protein = parts[0]
        region  = parts[1].replace("_", "-")
    elif (number_of_underscores == 2):
        parts   = protein_id.split('_')
        protein = parts[0]
        region  = f"{parts[1]}-{parts[2]}"
    elif (number_of_underscores == 3):
        parts        = protein_id.split('_')
        protein      = parts[0]
        start        = int(parts[1])
        region_parts = protein_id.split('/')[1].split('_')
        region       = f"{start + int(region_parts[0]) - 1}-{start + int(region_parts[1]) - 1}"

    return protein, region

def process_file(file_path, output_fasta):
    basename = os.path.basename(file_path).split('.')[0]

    alignment = AlignIO.read(file_path, "stockholm")
    first_record = alignment[0]
    sequence = ''.join([c.upper() for c in first_record.seq if c.isalpha()])
    sequence_id, region = parse_protein_region(first_record.id)

    with open(output_fasta, 'a') as out_f:
        if (region != "-"):
            out_f.write(f">{basename}-{sequence_id}_{region}\n{sequence}\n")
        else:
            out_f.write(f">{basename}-{sequence_id}\n{sequence}\n")

def main(input_msa_folder, output_fasta):
    for file in os.listdir(input_msa_folder):
        process_file(os.path.join(input_msa_folder, file), output_fasta)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_first_stockholm_sequences.py <input_msa_folder> <output_fasta>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
