import re
import sys
import bz2

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python3 bin/keep_unannotated_fasta.py test-data/sequence_explorer_protein_test.csv.bz2 data/output/unannotated.fasta")
    sys.exit(1)

input_filename = sys.argv[1]
output_filename = sys.argv[2]

# Compile the regex pattern outside of the loop for efficiency
pattern = re.compile(r'""p""')

with bz2.open(input_filename, 'rt', encoding='utf-8') as infile, open(output_filename, 'w') as outfile:
    # Skip the header line
    next(infile)
    
    for line in infile:
        # Split the line into columns
        mgyp, sequence, _, _, metadata_str = line.strip().split(',', 4)
        
        # Check for the "p" key using regex
        if not pattern.search(metadata_str):
            # Write to the output file in FASTA format
            outfile.write(f'>{mgyp}\n{sequence}\n')
