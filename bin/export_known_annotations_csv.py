from Bio import SeqIO
import csv
import argparse
import os

def process_record(record, writer, first_rep_id):
    split_id = record.id.split('|')
    if len(split_id) >= 3:
        rep_id, uniprot_name = split_id[1:3]
        split_description = record.description.split('OS=')
        if len(split_description[0].split(maxsplit=1)) > 1:
            uniprot_descr = split_description[0].split(maxsplit=1)[1].strip()
        else:
            uniprot_descr = None
        writer.writerow([first_rep_id, uniprot_name, uniprot_descr, "Uniprot SP", True])

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as fasta, open(output_file, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        records = SeqIO.parse(fasta, "fasta")
        first_record = next(records)
        first_rep_id = first_record.id.split('|')[1] if '|' in first_record.id else first_record.id
        if first_record.id.startswith('sp') and '|' in first_record.id:
            process_record(first_record, writer, first_rep_id)
        for record in records:
            if record.id.startswith('sp') and '|' in record.id:
                process_record(record, writer, first_rep_id)

    # check if file is empty and remove it
    if os.path.getsize(output_file) == 0:
        os.remove(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a FASTA file and write results to CSV.")
    parser.add_argument('input_file', type=str, help='The input FASTA file.')
    parser.add_argument('output_file', type=str, help='The output CSV file.')
    args = parser.parse_args()

    process_fasta(args.input_file, args.output_file)
