import pandas as pd
import sys

def find_duplicates(file_path, column_index):
    # Read the dataset into a pandas DataFrame
    data = pd.read_csv(file_path, delimiter='\t')
    
    # Get the column name based on the index
    column_name = data.columns[column_index]
    
    # Find the duplicate values in the specified column
    duplicates = data[data.duplicated(subset=column_name, keep=False)]
    
    return duplicates

if __name__ == '__main__':
    # Check if the required arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <file_path> <column_index>")
        sys.exit(1)
    
    # Extract the arguments
    file_path = sys.argv[1]
    column_index = int(sys.argv[2])
    
    # Call the function to find duplicates
    duplicates = find_duplicates(file_path, column_index)
    
    if duplicates.empty:
        print(f"No duplicates found in column at index {column_index}.")
    else:
        print(f"Duplicates found in column at index {column_index}:")
        print(duplicates.iloc[:, column_index].tolist())
