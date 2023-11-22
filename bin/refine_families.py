import sys

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 refine_families.py [Families TSV] [FASTA file] [Minimum number of family members] [Output file]")
        sys.exit(1)

    families_tsv = sys.argv[1]
    fasta_file = sys.argv[2]
    minimum_members = sys.argv[3]
    output_file = sys.argv[4]

    print("First message to console.")
    print("Second message to console.")

    with open("refined_families.tsv", "w") as file:
        file.write("First line of text.\n")

if __name__ == "__main__":
    main()
