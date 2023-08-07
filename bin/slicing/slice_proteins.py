import sys
from Bio import SeqIO
import pandas as pd

def slice_protein_sequence(fasta_file, annotation_file, min_slice_length):
    protein_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    annotations = pd.read_csv(annotation_file)
    annotations = annotations.sort_values(['ProteinID', 'Start'])

    # Create an empty list to store resulting sequences
    results = []
    previous_protein_id = None
    previous_end = 0

    # Loop over rows in the dataframe
    for index, row in annotations.iterrows():
        if previous_protein_id is not None and previous_protein_id != row['ProteinID']:
            # Check if there is any remaining sequence after the last annotation of the previous protein
            if len(protein_records[previous_protein_id]) - previous_end > min_slice_length:
                results.append((f"{previous_protein_id}_{previous_end + 1}_{len(protein_records[previous_protein_id])}",
                                str(protein_records[previous_protein_id].seq[previous_end:len(protein_records[previous_protein_id])])))
                previous_end = 0

        if row['Start'] - previous_end > min_slice_length:
            # Append a tuple with protein id and sequence to the results list
            results.append((f"{row['ProteinID']}_{previous_end + 1}_{row['Start'] - 1}",
                            str(protein_records[row['ProteinID']].seq[previous_end:row['Start'] - 1])))

        previous_end = row['End']
        previous_protein_id = row['ProteinID']

    # Check if there is any remaining sequence after the last annotation of the last protein
    if len(protein_records[previous_protein_id]) - previous_end > min_slice_length:
        results.append((f"{previous_protein_id}_{previous_end + 1}_{len(protein_records[previous_protein_id])}",
                        str(protein_records[previous_protein_id].seq[previous_end:len(protein_records[previous_protein_id])])))

    return results


def main():
    fasta_file = sys.argv[1]
    annotation_file = sys.argv[2]
    min_slice_length = int(sys.argv[3])
    output_file = sys.argv[4]

    slices = slice_protein_sequence(fasta_file, annotation_file, min_slice_length)

    with open(output_file, 'w') as f:
        for slice_info in slices:
            f.write(f">{slice_info[0]}\n{slice_info[1]}\n")
    
if __name__ == "__main__":
    main()

