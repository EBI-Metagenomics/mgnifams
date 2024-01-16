import sys
import os
import glob
import csv
from Bio import SeqIO

def initiate_output_csvs(mgnifams_out_dir, output_dir):
    mgnifam_headers          = ['id', 'size', 'protein_rep', 'mask', 'cif_file', 'msa_file', 'hmm_file']
    mgnifam_proteins_headers = ['mgnifam_id', 'protein', 'mask']
    mgnifam_pfams_headers    = ['mgnifam_id', 'rank', 'pfam_id', 'pfam_hit', 'query_hmm_range', 'template_hmm_range', 'e_value']
    mgnifam_folds_headers    = ['mgnifam_id', 'target_structure', 'aligned_length', 'query_start', 'query_end', 'target_start', 'target_end', 'e_value']

    mgnifam_csv_path          = os.path.join(output_dir, 'mgnifam.csv')
    mgnifam_proteins_csv_path = os.path.join(output_dir, 'mgnifam_proteins.csv')
    mgnifam_pfams_csv_path    = os.path.join(output_dir, 'mgnifam_pfams.csv')
    mgnifam_folds_csv_path    = os.path.join(output_dir, 'mgnifam_folds.csv')

    with open(mgnifam_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(mgnifam_headers)
    with open(mgnifam_proteins_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(mgnifam_proteins_headers)
    with open(mgnifam_pfams_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(mgnifam_pfams_headers)
    with open(mgnifam_folds_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(mgnifam_folds_headers)

def get_number_of_families(mgnifams_out_dir):
    msa_directory = os.path.join(mgnifams_out_dir, 'families', 'msa')
    msa_files = glob.glob(os.path.join(msa_directory, '*'))
    num_mgnifams = len(msa_files)
    return num_mgnifams

def parse_cif(mgnifam_id, mgnifams_out_dir):
    cif_directory = os.path.join(mgnifams_out_dir, 'cif')
    cif_files = glob.glob(os.path.join(cif_directory, '**', mgnifam_id + '_*'), recursive=True)
    cif_filename = os.path.basename(cif_files[0]) if cif_files else None

    filename_no_ext = cif_filename.split('.')[0]
    first_split = filename_no_ext.split('-')
    first_split_first_part = first_split[0]
    first_split_second_part = first_split[1]
    size = first_split_first_part.split('_')[1]
    protein_parts = first_split_second_part.split('_')
    protein_rep = protein_parts[0]
    if (len(protein_parts) == 1):
        mask = "-"
    elif (len(protein_parts) == 3):
        mask = f"{protein_parts[1]}-{protein_parts[2]}"
    elif (len(protein_parts) == 5):
        mask = f"{int(protein_parts[1]) + int(protein_parts[3]) - 1}-{int(protein_parts[1]) + int(protein_parts[4]) - 1}"

    return size, protein_rep, mask, cif_filename

def write_mgnifam(i, mgnifams_out_dir, output_dir):
    mgnifam_id = f"mgnfam{i}"
    size, protein_rep, mask, cif_file = parse_cif(mgnifam_id, mgnifams_out_dir)
    msa_file = f"{mgnifam_id}_{size}.fa"
    hmm_file = f"{mgnifam_id}_{size}.hmm"
    
    mgnifam_csv_path = os.path.join(output_dir, 'mgnifam.csv')
    with open(mgnifam_csv_path, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([mgnifam_id, size, protein_rep, mask, cif_file, msa_file, hmm_file]) 

def parse_protein_mask(protein_id):
    number_of_underscores = protein_id.count('_')
    if (number_of_underscores == 0):
        protein = protein_id
        mask    = "-"
    elif (number_of_underscores == 1):
        parts   = protein_id.split('/')
        protein = parts[0]
        mask    = parts[1].replace("_", "-")
    elif (number_of_underscores == 2):
        parts   = protein_id.split('_')
        protein = parts[0]
        mask    = f"{parts[1]}-{parts[2]}"
    elif (number_of_underscores == 3):
        parts      = protein_id.split('_')
        protein    = parts[0]
        start      = int(parts[1])
        mask_parts = protein_id.split('/')[1].split('_')
        mask       = f"{start + int(mask_parts[0]) - 1}-{start + int(mask_parts[1]) - 1}"

    return protein, mask

def write_mgnifam_proteins(mgnifams_out_dir, output_dir):
    refined_families_file = os.path.join(mgnifams_out_dir, 'families', 'updated_refined_families.tsv')
    mgnifam_proteins_csv_path = os.path.join(output_dir, 'mgnifam_proteins.csv')

    with open(refined_families_file, 'r') as file, open(mgnifam_proteins_csv_path, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file)

        for line in file:
            mgnifam_id, protein_mask = line.strip().split('\t')
            protein, mask = parse_protein_mask(protein_mask)
            writer.writerow([mgnifam_id, protein, mask])

def write_mgnifam_pfams(mgnifams_out_dir, output_dir):
    pfam_hits_directory = os.path.join(mgnifams_out_dir, 'hh', 'hits')
    hits_files = glob.glob(os.path.join(pfam_hits_directory, '*'))

    mgnifam_pfams_csv_path = os.path.join(output_dir, 'mgnifam_pfams.csv')

    for file_path in hits_files:
        with open(file_path, 'r') as file:
            hits_data = []
            mgnifam_id = os.path.basename(file_path).split('_')[0]

            for line in file:
                name = line[4:34].strip()
                pfam_id = name.split(';')[0].strip().split('.')[0]
                
                hit = {
                    'mgnifam_id': mgnifam_id,
                    'rank': line[0:3].strip(),
                    'pfam_id': pfam_id,
                    'name': name,
                    'query_hmm': line[75:83].strip(),
                    'template_hmm': line[84:99].strip(),
                    'e_value': line[41:48].strip()
                }
                hits_data.append(hit)

            with open(mgnifam_pfams_csv_path, 'a', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=hit.keys())
                if os.path.getsize(mgnifam_pfams_csv_path) == 0:  # If file is empty, write header
                    writer.writeheader()
                writer.writerows(hits_data)

def write_mgnifam_folds(mgnifams_out_dir, output_dir):
    foldseek_directory = os.path.join(mgnifams_out_dir, 'foldseek')
    foldseek_files = glob.glob(os.path.join(foldseek_directory, '*'))

    mgnifam_folds_csv_path = os.path.join(output_dir, 'mgnifam_folds.csv')

    for filepath in foldseek_files:
        filename = os.path.basename(filepath)
        if filename.startswith(('alphafold_', 'esm_', 'pdb_')):
            with open(filepath, 'r') as f:
                structural_annotations = []
                for file_line in f:
                    parts = file_line.strip().split('\t')
                    mgnifam_id = parts[0].split('-')[0].split('_')[0]
                    target_structure_identifier = parts[1]

                    annotation = {
                        'mgnifam_id': mgnifam_id,
                        'target_structure_identifier': target_structure_identifier,
                        'aligned_length': int(parts[3]),
                        'query_start': int(parts[6]),
                        'query_end': int(parts[7]),
                        'target_start': int(parts[8]),
                        'target_end': int(parts[9]),
                        'e_value': float(parts[10])
                    }
                    structural_annotations.append(annotation)

            # Write annotations to CSV file
            with open(mgnifam_folds_csv_path, 'a', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=annotation.keys())
                if os.path.getsize(mgnifam_folds_csv_path) == 0:  # If file is empty, write header
                    writer.writeheader()
                writer.writerows(structural_annotations)

def main():
    if len(sys.argv) != 2:
        print("Usage: python export_mgnifams_csv.py <mgnifams_out_dir>")
        sys.exit(1)

    mgnifams_out_dir = sys.argv[1]
    output_dir       = os.path.join(mgnifams_out_dir, "tables")
    os.makedirs(output_dir, exist_ok=True)

    initiate_output_csvs(mgnifams_out_dir, output_dir)
    num_mgnifams = get_number_of_families(mgnifams_out_dir)
    for i in range(1, num_mgnifams + 1):
        write_mgnifam(i, mgnifams_out_dir, output_dir)
        
    write_mgnifam_proteins(mgnifams_out_dir, output_dir)
    write_mgnifam_pfams(mgnifams_out_dir, output_dir)
    write_mgnifam_folds(mgnifams_out_dir, output_dir)

if __name__ == "__main__":
    main()
