#!/usr/bin/env python3

import argparse
import os

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
        "--output_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output file with pre-redundant families TSV edgelist.",
    )
    return parser.parse_args(args)
    
def pool_directory(input_dir, output_filename, splitChar):
    
    with open(output_filename, 'w') as outfile:
        for filename in sorted(os.listdir(input_dir)):
            filepath = os.path.join(input_dir, filename)
            with open(filepath, 'r') as infile:
                base_filename = os.path.splitext(filename)[0]
                for line in infile:
                    value = line.strip()
                    if value:  # Check if the line is not empty
                        if splitChar != "":
                            first_element       = value.split(splitChar)[0]
                            fam_id              = base_filename + '_' + first_element
                            remaining_elements  = splitChar.join(value.split(splitChar)[1:])
                            concatenated_string = fam_id + splitChar + remaining_elements

                            outfile.write(f"{concatenated_string}\n")

def main(args=None):
    args = parse_args(args)

    global arg_families_dir, arg_out_file

    arg_families_dir = args.input_dir
    arg_out_file = args.output_file

    pool_directory(arg_families_dir, arg_out_file, "\t")

if __name__ == "__main__":
    main()