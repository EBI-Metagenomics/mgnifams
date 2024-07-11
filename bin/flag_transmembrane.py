#!/usr/bin/env python3

import sys
import re

def parse_args():
    global arg_fraction, arg_gff3, arg_outfile
        
    if not (len(sys.argv) == 4):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_fraction = float(sys.argv[1])
    arg_gff3     = sys.argv[2]
    arg_outfile  = sys.argv[3]

def flag_transmembrane_proteins():
    results = {}

    with open(arg_gff3, 'r') as file:
        lines = file.readlines()

    # Regular expression to match headers
    header_re = re.compile(r'# (\S+) Length: (\d+)')

    current_id = None
    total_length = 0
    tmhelix_length = 0

    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            match = header_re.match(line)
            if match:
                if current_id:
                    results[current_id] = tmhelix_length / total_length if total_length > 0 else 0
                current_id = match.group(1)
                total_length = int(match.group(2))
                tmhelix_length = 0
        elif line.startswith('//'):
            if current_id:
                results[current_id] = tmhelix_length / total_length if total_length > 0 else 0
            current_id = None
        else:
            parts = line.split('\t')
            if parts[1] == 'TMhelix' or parts[1] == 'Beta sheet':
                tmhelix_length += int(parts[3]) - int(parts[2]) + 1

    # If the last ID was not processed yet
    if current_id:
        results[current_id] = tmhelix_length / total_length if total_length > 0 else 0

    with open(arg_outfile, 'w') as f:
        for id, fraction in results.items():
            print(f"{id}: {fraction:.4f}")
            if fraction >= arg_fraction:
                f.write(f"{id}\n")

def main():
    parse_args()
    flag_transmembrane_proteins()

if __name__ == "__main__":
    main()
