import pandas as pd
import argparse
import os

def create_csv(input_tsv, output_csv):
    try:
        # Read the TSV file into a pandas DataFrame
        df = pd.read_csv(input_tsv, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        # Create an empty DataFrame if there's an error reading the input file
        df = pd.DataFrame()

    if not df.empty:
        # Concatenate columns 5, 11, and 12
        df['concat'] = df.iloc[:,5].astype(str) + '_' + df.iloc[:,11].astype(str) + '_' + df.iloc[:,12].astype(str)

        # Select the required columns and rename them
        df_out = df[[0,4,'concat',3]].copy()
        df_out['isKnown'] = False
        df_out.columns = ['FamilyID', 'Annotation', 'Description', 'Source', 'IsKnown']

        # Write the DataFrame to a CSV file
        df_out.to_csv(output_csv, index=False, header=False)
    else:
        # Create an empty file if the DataFrame is empty
        open(output_csv, 'w').close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a CSV file with specific columns from a TSV file')
    parser.add_argument('input_tsv', type=str, help='Path to the input TSV file')
    parser.add_argument('output_csv', type=str, help='Path to the output CSV file')

    args = parser.parse_args()

    create_csv(args.input_tsv, args.output_csv)