#!/usr/bin/env python3

import sys
import os

def parse_args():
    global arg_families_dir, arg_outfile
        
    if not (len(sys.argv) == 3):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_families_dir = sys.argv[1]
    arg_outfile      = sys.argv[2]
    
def pool_directory(input_dir, output_file, splitChar):
    
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(arg_families_dir):
            filepath = os.path.join(arg_families_dir, filename)
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
                        else:
                            pass

                        outfile.write(f"{concatenated_string}\n")

def main():
    parse_args()    
    pool_directory("refined_families", arg_outfile, "\t")

if __name__ == "__main__":
    main()
