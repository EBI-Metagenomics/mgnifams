import sys

def read_ids(id_file):
    with open(id_file, 'r') as file:
        return set(line.strip() for line in file)

def find_fasta_by_id(id_set, fasta_file, output_file):
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        write_entry = False
        for line in fasta:
            if line.startswith('>'):
                entry_id = line.split()[0][1:]  # Extract the ID part
                write_entry = entry_id in id_set
            if write_entry:
                output.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 find_fasta_by_id.py [ID file] [FASTA file] [Output file]")
        sys.exit(1)

    id_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    id_set = read_ids(id_file)
    find_fasta_by_id(id_set, fasta_file, output_file)

if __name__ == "__main__":
    main()