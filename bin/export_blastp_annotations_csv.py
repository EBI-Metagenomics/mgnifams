import pandas as pd
import argparse
import os

def create_csv(input_tsv, output_csv):
    try:
        df = pd.read_csv(input_tsv, delimiter='\t', header=None)
    except pd.errors.EmptyDataError:
        # Create an empty DataFrame if there's an error reading the input file
        df = pd.DataFrame()

    if not df.empty:
        # Split third column using "|" as separator and keep the second part
        try:
            df[1] = df[1].apply(lambda x: x.split('|')[2])
        except IndexError:
            print(f"Error: Unable to split column values as expected in {input_tsv}")
            return
        
        # Create new columns with static values
        df['None'] = None
        df['Uniprot SP'] = 'Uniprot SP'
        df['False'] = False

        # Select the required columns and rename them
        df_out = df[[0, 1, 'None', 'Uniprot SP', 'False']].copy()
        df_out.columns = ['FamilyID', 'Annotation', 'Description', 'Source', 'IsKnown']

        df_out.to_csv(output_csv, sep=',', index=False, header=False)
    else:
        # Create an empty file if the DataFrame is empty
        open(output_csv, 'w').close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a CSV file with specific columns from a TSV file')
    parser.add_argument('input_tsv', type=str, help='Path to the input TSV file')
    parser.add_argument('output_csv', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    create_csv(args.input_tsv, args.output_csv)