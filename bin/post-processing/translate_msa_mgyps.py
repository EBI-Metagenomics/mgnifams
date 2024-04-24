import sys
import os

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

def format_protein_name(raw_name):
    """
    Formats the protein name by appending zeros in front to make it 12 characters,
    and then adds 'MGYP' as a prefix.
    """
    formatted_name = raw_name.zfill(12)  # Append zeros to make it 12 characters
    return "MGYP" + formatted_name

def process_file(input_file, output_file):
    print(input_file)
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                protein_id = line.strip()[1:]  # Remove the '>'
                new_protein, new_region = parse_protein_region(protein_id)
                formatted_name = format_protein_name(new_protein)
                f_out.write(f'>{formatted_name}/{new_region}\n')
            else:
                f_out.write(line)

    os.replace(output_file, input_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python bin/post-processing/translate_msa_mgyps.py input_folder") # msa or seed_msa paths
    else:
        input_folder = sys.argv[1]

        for file_name in os.listdir(input_folder):
            if file_name.endswith(".fas"):
                input_file = os.path.join(input_folder, file_name)
                output_file = "temp_out.fas"
                process_file(input_file, output_file)
