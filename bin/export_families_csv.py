import argparse
import csv

def export_families_csv(reps_file, csv_file, mode):
    with open(reps_file, 'r') as reps:
        with open(csv_file, 'w', newline='') as csvfile:
            fieldnames = ['ID', 'IsKnown', 'MSA', 'HMM']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for rep in reps:
                rep = rep.strip()  # remove trailing newline
                writer.writerow({
                    'ID': rep,
                    'IsKnown': 1 if mode == 'known' else 0,
                    'MSA': f'data/output/mafft/{rep}_msa.fa',
                    'HMM': f'data/output/hmm/build/{rep}_family.hmm'
                })

def main():
    parser = argparse.ArgumentParser(description='Export MGnifams DB Families CSV file for import.')
    parser.add_argument('reps_file', help='The path to the file containing the cluster rep names.')
    parser.add_argument('csv_file', help='The path to the output CSV file.')
    parser.add_argument('mode', help='known or unknown')
    args = parser.parse_args()

    export_families_csv(args.reps_file, args.csv_file, args.mode)

if __name__ == "__main__":
    main()
