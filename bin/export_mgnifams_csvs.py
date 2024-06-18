import sys
import os
import glob
import csv
import pandas as pd

def initiate_output_csvs(mgnifams_out_dir, output_dir):
    mgnifam_headers = ['id', 'family_size', 'protein_rep', 'rep_region', 'rep_length', 'plddt', 'ptm', 'converged', \
        'cif_file', 'seed_msa_file', 'msa_file', 'hmm_file', 'rf_file', 'biomes_file', 'domain_architecture_file']
    mgnifam_proteins_headers = ['mgnifam_id', 'protein', 'region']
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

def read_family_metadata(mgnifams_out_dir, families_dir_name):
    column_names = ['family', 'size', 'protein_rep', 'region']
    family_metadata_df = pd.read_csv(os.path.join(mgnifams_out_dir, families_dir_name, "family_metadata.csv"), header=None, names=column_names)
    return family_metadata_df

def read_structure_scores(mgnifams_out_dir):
    column_names = ['family', 'rep_length', 'plddt', 'ptm']
    structure_scores_df = pd.read_csv(os.path.join(mgnifams_out_dir, "structures/pdb_scores.csv"), header=None, names=column_names)
    return structure_scores_df

def get_converged_families(mgnifams_out_dir, families_dir_name):
    with open(os.path.join(mgnifams_out_dir, families_dir_name, 'converged_families.txt'), 'r') as file:
        converged_families = {line.strip() for line in file}

        return converged_families

def is_converged(fam, converged_families):
    return str(fam) in converged_families

def write_mgnifam(family_id, mgnifams_out_dir, output_dir, family_metadata_df, structure_scores_df, converged_families):
    family_size, protein_rep, region = family_metadata_df[family_metadata_df['family'] == family_id][['size', 'protein_rep', 'region']].values[0]
    rep_length, plddt, ptm = structure_scores_df[structure_scores_df['family'] == family_id][['rep_length', 'plddt', 'ptm']].values[0]
    converged = is_converged(family_id, converged_families)
    
    cif_file                 = f"{family_id}.cif"
    seed_msa_file            = f"{family_id}.fas"
    msa_file                 = f"{family_id}.fas"
    hmm_file                 = f"{family_id}.hmm"
    rf_file                  = f"{family_id}.txt"
    biomes_file              = f"{family_id}.csv"
    domain_architecture_file = f"{family_id}.json"
    
    mgnifam_csv_path = os.path.join(output_dir, 'mgnifam.csv')
    with open(mgnifam_csv_path, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([family_id, family_size, protein_rep, region, rep_length, plddt, ptm, converged,
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

def write_mgnifam_proteins(mgnifams_out_dir, families_dir_name, output_dir):
    refined_families_file = os.path.join(mgnifams_out_dir, families_dir_name, 'refined_families.tsv')
    mgnifam_proteins_csv_path = os.path.join(output_dir, 'mgnifam_proteins.csv')

    with open(refined_families_file, 'r') as file, open(mgnifam_proteins_csv_path, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file)

        for line in file:
            mgnifam_id, protein_region = line.strip().split('\t')
            mgnifam_id = mgnifam_id.replace("mgnifam", "")
            protein, region = parse_protein_region(protein_region)
            writer.writerow([mgnifam_id, protein, region])

def write_mgnifam_pfams(mgnifams_out_dir, output_dir):
    pfam_hits_directory = os.path.join(mgnifams_out_dir, 'hh', 'pfam_hits')
    hits_files = glob.glob(os.path.join(pfam_hits_directory, '*'))

    mgnifam_pfams_csv_path = os.path.join(output_dir, 'mgnifam_pfams.csv')

    for file_path in hits_files:
        with open(file_path, 'r') as file:
            hits_data = []
            mgnifam_id = os.path.basename(file_path).split('.')[0]

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
    foldseek_hits_file = os.path.join(mgnifams_out_dir, 'structures/foldseek/foldseek_hits.tsv')
    mgnifam_folds_csv_path = os.path.join(output_dir, 'mgnifam_folds.csv')

    annotation = {}
    with open(foldseek_hits_file, 'r') as f:
        structural_annotations = []
        for file_line in f:
            parts = file_line.strip().split('\t')
            mgnifam_id = parts[0].split('.')[0]
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
    if len(sys.argv) != 3:
        print("Usage: python bin/export_mgnifams_csvs.py <mgnifams_out_dir> <families_dir_name>")
        sys.exit(1)

    mgnifams_out_dir  = sys.argv[1]
    families_dir_name = sys.argv[2]
    output_dir        = 'tables'
    os.makedirs(output_dir, exist_ok=True)

    initiate_output_csvs(mgnifams_out_dir, output_dir)
    family_metadata_df  = read_family_metadata(mgnifams_out_dir, families_dir_name)
    structure_scores_df = read_structure_scores(mgnifams_out_dir)
    converged_families  = get_converged_families(mgnifams_out_dir, families_dir_name)

    for family_id in family_metadata_df['family']:
        write_mgnifam(family_id, mgnifams_out_dir, output_dir, family_metadata_df, structure_scores_df, converged_families)
        
    write_mgnifam_proteins(mgnifams_out_dir, families_dir_name, output_dir)
    write_mgnifam_pfams(mgnifams_out_dir, output_dir)
    write_mgnifam_folds(mgnifams_out_dir, output_dir)

if __name__ == "__main__":
    main()
