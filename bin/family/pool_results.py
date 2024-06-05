import sys
import os
import shutil
import json

def parse_args():
    global arg_families_dir
        
    if not (len(sys.argv) == 2):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_families_dir = sys.argv[1]

def create_mapping_dict():
    filenames = os.listdir(os.path.join(arg_families_dir, 'rf'))
    filenames_without_extension = [os.path.splitext(filename)[0] for filename in filenames]
    filenames_without_extension.sort()
    family_to_id = {filename: idx + 1 for idx, filename in enumerate(filenames_without_extension)}

    return family_to_id
    
def pool_directory(input_dir, output_filename, splitChar):
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file    = os.path.join(outdir, output_filename)    
    
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(path_to_folder):
            filepath = os.path.join(path_to_folder, filename)
            with open(filepath, 'r') as infile:
                base_filename = os.path.splitext(filename)[0]
                for line in infile:
                    value = line.strip()
                    if value:  # Check if the line is not empty
                        if splitChar != "":
                            first_element       = value.split(splitChar)[0]
                            remaining_elements  = splitChar.join(value.split(splitChar)[1:])
                            family_id           = str(family_to_id[base_filename + '_' + first_element])
                            concatenated_string = family_id + splitChar + remaining_elements
                        else: # converged_families file, only one value element
                            concatenated_string = str(family_to_id[base_filename + '_' + value])

                        outfile.write(f"{concatenated_string}\n")

def pool_clusters_directory(input_dir, output_filename):
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file    = os.path.join(outdir, output_filename)
    
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(path_to_folder):
            filepath = os.path.join(path_to_folder, filename)
            with open(filepath, 'r') as infile:
                content = infile.read().strip()  # Read and strip leading/trailing whitespace
                if content:  # Only write if the file is not empty
                    outfile.write(content + '\n')

def translate_directory(input_dir):
    input_folder  = os.path.join(arg_families_dir, input_dir)
    os.makedirs(os.path.join(outdir, input_dir), exist_ok=True)
    output_folder = os.path.join(outdir, input_dir)

    filenames = os.listdir(input_folder)
    for filename in filenames:
        basename         = os.path.splitext(filename)[0]
        source_path      = os.path.join(input_folder, filename)
        new_filename     = f'{family_to_id[basename]}{os.path.splitext(filename)[1]}'
        destination_path = os.path.join(output_folder, new_filename)
        
        shutil.copy(source_path, destination_path)
        
def main():
    parse_args()

    global outdir, family_to_id

    outdir = "families_pooled"
    os.makedirs(outdir, exist_ok=True)

    family_to_id = create_mapping_dict()

    pool_directory("family_metadata", "family_metadata.csv", ",")
    pool_directory("refined_families", "refined_families.tsv", "\t")
    pool_directory("converged_families", "converged_families.txt", "")
    pool_clusters_directory("successful_clusters", "successful_clusters.txt")
    pool_clusters_directory("discarded_clusters", "discarded_clusters.txt")

    translate_directory('rf')
    translate_directory('hmm')
    translate_directory('msa_sto')
    translate_directory('seed_msa_sto')
    translate_directory('domtblout')
    
    json_mapping = 'family_to_id.json'
    output_file  = os.path.join(outdir, json_mapping)
    with open(output_file, 'w') as f:
        json.dump(family_to_id, f)

if __name__ == "__main__":
    main()