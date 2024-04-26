import sys
import os
import glob
import csv

def initiate_output_csvs(khalifams_out_dir, output_dir):
    khalifam_headers          = ['id', 'family_size', 'protein_rep', 'rep_region', 'converged', 'cif_file', 'seed_msa_file', 'msa_file', 'hmm_file', 'rf_file', 'biomes_file', 'domain_architecture_file']
    khalifam_proteins_headers = ['khalifam_id', 'protein', 'region']
    khalifam_pfams_headers    = ['khalifam_id', 'rank', 'pfam_id', 'pfam_hit', 'query_hmm_range', 'template_hmm_range', 'e_value']
    khalifam_folds_headers    = ['khalifam_id', 'target_structure', 'aligned_length', 'query_start', 'query_end', 'target_start', 'target_end', 'e_value']

    khalifam_csv_path          = os.path.join(output_dir, 'khalifam.csv')
    khalifam_proteins_csv_path = os.path.join(output_dir, 'khalifam_proteins.csv')
    khalifam_pfams_csv_path    = os.path.join(output_dir, 'khalifam_pfams.csv')
    khalifam_folds_csv_path    = os.path.join(output_dir, 'khalifam_folds.csv')

    with open(khalifam_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(khalifam_headers)
    with open(khalifam_proteins_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(khalifam_proteins_headers)
    with open(khalifam_pfams_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(khalifam_pfams_headers)
    with open(khalifam_folds_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(khalifam_folds_headers)

def get_number_of_families(khalifams_out_dir):
    msa_directory = os.path.join(khalifams_out_dir, 'families', 'seed_msa')
    msa_files = glob.glob(os.path.join(msa_directory, '*'))
    num_khalifams = len(msa_files)
    return num_khalifams

def get_converged_families(khalifams_out_dir):
    with open(os.path.join(khalifams_out_dir, 'families', 'updated_converged_families.txt'), 'r') as file:
        converged_families = {line.strip() for line in file}

        return converged_families

def parse_cif(khalifam_id, khalifams_out_dir):
    cif_directory = os.path.join(khalifams_out_dir, 'cif')
    cif_files = glob.glob(os.path.join(cif_directory, '**', khalifam_id + '_*'), recursive=True)
    cif_filename = os.path.basename(cif_files[0]) if cif_files else None

    filename_no_ext = cif_filename.split('.')[0]
    first_split = filename_no_ext.split('-')
    first_split_first_part = first_split[0]
    first_split_second_part = first_split[1]
    family_size = first_split_first_part.split('_')[1]
    protein_parts = first_split_second_part.split('_')
    protein_rep = protein_parts[0]
    if (len(protein_parts) == 1):
        region = "-"
    elif (len(protein_parts) == 3):
        region = f"{protein_parts[1]}-{protein_parts[2]}"
    elif (len(protein_parts) == 5):
        region = f"{int(protein_parts[1]) + int(protein_parts[3]) - 1}-{int(protein_parts[1]) + int(protein_parts[4]) - 1}"

    return family_size, protein_rep, region, cif_filename

def is_converged(fam, converged_families):
    return fam in converged_families

def write_khalifam(i, khalifams_out_dir, output_dir, converged_families):
    khalifam_id = f"khalifam{i}"
    family_size, protein_rep, region, cif_file = parse_cif(khalifam_id, khalifams_out_dir)
    converged = is_converged(f"khalifam{i}", converged_families)
    seed_msa_file = f"{khalifam_id}_{family_size}.fas"
    msa_file = seed_msa_file
    hmm_file = f"{khalifam_id}_{family_size}.hmm"
    rf_file = f"{khalifam_id}_{family_size}.txt"
    biomes_file = f"{khalifam_id}_b_counts.csv"
    domain_architecture_file = f"{khalifam_id}_domains.json"
    
    khalifam_csv_path = os.path.join(output_dir, 'khalifam.csv')
    with open(khalifam_csv_path, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([i, family_size, protein_rep, region, converged,
            cif_file, seed_msa_file, msa_file, hmm_file, rf_file, biomes_file, domain_architecture_file]) 

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

def write_khalifam_proteins(khalifams_out_dir, output_dir):
    refined_families_file = os.path.join(khalifams_out_dir, 'families', 'updated_refined_families.tsv')
    khalifam_proteins_csv_path = os.path.join(output_dir, 'khalifam_proteins.csv')

    with open(refined_families_file, 'r') as file, open(khalifam_proteins_csv_path, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file)

        for line in file:
            khalifam_id, protein_region = line.strip().split('\t')
            khalifam_id = khalifam_id.replace("khalifam", "")
            protein, region = parse_protein_region(protein_region)
            writer.writerow([khalifam_id, protein, region])

def write_khalifam_pfams(khalifams_out_dir, output_dir):
    pfam_hits_directory = os.path.join(khalifams_out_dir, 'hh', 'hits')
    hits_files = glob.glob(os.path.join(pfam_hits_directory, '*'))

    khalifam_pfams_csv_path = os.path.join(output_dir, 'khalifam_pfams.csv')

    for file_path in hits_files:
        with open(file_path, 'r') as file:
            hits_data = []
            khalifam_id = os.path.basename(file_path).split('_')[0]
            khalifam_id = khalifam_id.replace("khalifam", "")

            for line in file:
                name = line[4:34].strip()
                pfam_id = name.split(';')[0].strip().split('.')[0]
                
                hit = {
                    'khalifam_id': khalifam_id,
                    'rank': line[0:3].strip(),
                    'pfam_id': pfam_id,
                    'name': name,
                    'query_hmm': line[75:83].strip(),
                    'template_hmm': line[84:99].strip(),
                    'e_value': line[41:48].strip()
                }
                hits_data.append(hit)

            with open(khalifam_pfams_csv_path, 'a', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=hit.keys())
                if os.path.getsize(khalifam_pfams_csv_path) == 0:  # If file is empty, write header
                    writer.writeheader()
                writer.writerows(hits_data)

def write_khalifam_folds(khalifams_out_dir, output_dir):
    foldseek_directory = os.path.join(khalifams_out_dir, 'foldseek')
    foldseek_files = glob.glob(os.path.join(foldseek_directory, '*'))
    khalifam_folds_csv_path = os.path.join(output_dir, 'khalifam_folds.csv')

    annotation = {}
    for filepath in foldseek_files:
        filename = os.path.basename(filepath)
        if filename.startswith(('alphafold_', 'esm_', 'pdb_')):
            with open(filepath, 'r') as f:
                structural_annotations = []
                for file_line in f:
                    parts = file_line.strip().split('\t')
                    khalifam_id = parts[0].split('-')[0].split('_')[0]
                    khalifam_id = khalifam_id.replace("khalifam", "")
                    target_structure_identifier = parts[1]

                    annotation = {
                        'khalifam_id': khalifam_id,
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
            with open(khalifam_folds_csv_path, 'a', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=annotation.keys())
                if os.path.getsize(khalifam_folds_csv_path) == 0:  # If file is empty, write header
                    writer.writeheader()
                writer.writerows(structural_annotations)

def main():
    if len(sys.argv) != 2:
        print("Usage: python bin/export_khalifams_csvs.py /home/vangelis/Desktop/Projects/khalifams-site-data_backup")
        sys.exit(1)

    khalifams_out_dir = sys.argv[1]
    output_dir        = "tables"
    os.makedirs(output_dir, exist_ok=True)

    initiate_output_csvs(khalifams_out_dir, output_dir)
    num_khalifams = get_number_of_families(khalifams_out_dir)
    converged_families = get_converged_families(khalifams_out_dir)
    
    for i in range(1, num_khalifams + 1):
        write_khalifam(i, khalifams_out_dir, output_dir, converged_families)
        
    write_khalifam_proteins(khalifams_out_dir, output_dir)
    write_khalifam_pfams(khalifams_out_dir, output_dir)
    write_khalifam_folds(khalifams_out_dir, output_dir)

if __name__ == "__main__":
    main()
