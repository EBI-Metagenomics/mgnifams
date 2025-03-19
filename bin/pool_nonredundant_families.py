#!/usr/bin/env python3

import sys
import argparse
import os
import shutil
import json

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--input_dir",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Folder of all families, including redundant ones.",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder non-redundant families.",
    )
    parser.add_argument(
        "-r",
        "--redundant",
        required=True,
        metavar="FILE",
        type=str,
        help="File with redundant family ids.",
    )
    parser.add_argument(
        "-i",
        "--iter",
        required=True,
        metavar="INT",
        type=int,
        help="Starting family id.",
    )
    return parser.parse_args(args)

def read_non_redundant_fam_ids(file_path):
    with open(file_path, 'r') as file:
        lines = file.read().splitlines()
    return lines

def create_mapping_dict():
    hmm_folder  = os.path.join(arg_families_dir, 'hmm')
    hmm_files = sorted([name for name in os.listdir(hmm_folder) if os.path.basename(name).split('.')[0] not in redundant_fam_ids])

    family_to_id = {os.path.basename(filename).split('.')[0]: idx + arg_starting_id for idx, filename in enumerate(hmm_files)}

    return family_to_id
    
def pool_directory(input_dir, output_filename, splitChar):
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file    = os.path.join(arg_out_dir, output_filename)    
    
    with open(output_file, 'w') as outfile:
        for filename in sorted(os.listdir(path_to_folder)):
            filepath = os.path.join(path_to_folder, filename)
            with open(filepath, 'r') as infile:
                base_filename = os.path.splitext(filename)[0]
                for line in infile:
                    value = line.strip()
                    if value:  # Check if the line is not empty
                        if splitChar != "":
                            first_element       = value.split(splitChar)[0]
                            fam_id              = base_filename + '_' + first_element
                            if (fam_id in redundant_fam_ids):
                                continue
                            remaining_elements  = splitChar.join(value.split(splitChar)[1:])
                            family_id           = str(family_to_id[fam_id])
                            concatenated_string = family_id + splitChar + remaining_elements
                        else: # converged_families file, only one value element
                            fam_id              = base_filename + '_' + value
                            if (fam_id in redundant_fam_ids):
                                continue
                            concatenated_string = str(family_to_id[fam_id])

                        outfile.write(f"{concatenated_string}\n")

def pool_clusters_directory(input_dir, output_filename):
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file    = os.path.join(arg_out_dir, output_filename)
    
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(path_to_folder):
            filepath = os.path.join(path_to_folder, filename)
            with open(filepath, 'r') as infile:
                content = infile.read().strip()  # Read and strip leading/trailing whitespace
                if content:  # Only write if the file is not empty
                    outfile.write(content + '\n')

def translate_directory(input_dir):
    input_folder  = os.path.join(arg_families_dir, input_dir)
    os.makedirs(os.path.join(arg_out_dir, input_dir), exist_ok=True)
    output_folder = os.path.join(arg_out_dir, input_dir)

    filenames = os.listdir(input_folder)
    for filename in filenames:
        basename         = os.path.splitext(filename)[0]
        if (basename in redundant_fam_ids):
            continue
        source_path      = os.path.join(input_folder, filename)
        new_filename     = f'{family_to_id[basename]}{os.path.splitext(filename)[1]}'
        destination_path = os.path.join(output_folder, new_filename)
        
        shutil.copy(source_path, destination_path)

# def translate_edgelist(file_path, out_path):
#     if os.path.getsize(file_path) == 0:
#         shutil.copy(file_path, out_path)
#         return

#     df = pd.read_csv(file_path, header=None)
#     df.columns = ['Col1', 'Col2', 'Col3']
#     df['Col1'] = df['Col1'].map(family_to_id)
#     df['Col2'] = df['Col2'].map(family_to_id)
#     df.to_csv(out_path, index=False, header=False)

def main(args=None):
    args = parse_args(args)

    global arg_families_dir, arg_out_dir, \
        arg_non_redundant_fam_ids_file, arg_starting_id, \
        redundant_fam_ids, family_to_id

    arg_families_dir = args.input_dir
    arg_out_dir = args.output_dir
    arg_non_redundant_fam_ids_file = args.redundant
    arg_starting_id = args.iter

    os.makedirs(arg_out_dir, exist_ok=True)

    redundant_fam_ids = read_non_redundant_fam_ids(arg_non_redundant_fam_ids_file)
    family_to_id      = create_mapping_dict()

    pool_directory("family_metadata", "family_metadata.csv", ",")
    # pool_directory("family_reps", "family_reps.fasta", "")
    # pool_directory("refined_families", "refined_families.tsv", "\t")
    # pool_directory("converged_families", "converged_families.txt", "")
    # pool_clusters_directory("successful_clusters", "successful_clusters.txt")
    # pool_clusters_directory("discarded_clusters", "discarded_clusters.txt")

    # translate_directory('rf')
    # translate_directory('hmm')
    # translate_directory('msa_sto')
    # translate_directory('seed_msa_sto')
    # translate_directory('domtblout')
    # # translate_edgelist(arg_similarity_edgelist, \
    # #     os.path.join(arg_out_dir, 'similarity_edgelist.csv'))
    
    # json_mapping = 'family_to_id.json'
    # output_file  = os.path.join(arg_out_dir, json_mapping)
    # with open(output_file, 'w') as f:
    #     json.dump(family_to_id, f)

if __name__ == "__main__":
    main()