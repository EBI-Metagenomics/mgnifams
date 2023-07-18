import pandas as pd
import sys

def count_unique_elements(file_path, column_index):
    # Read the dataset into a pandas DataFrame
    data = pd.read_csv(file_path, delimiter="\t")
    
    # Get the column name based on the index
    column_name = data.columns[column_index]
    
    # Count the unique elements in the specified column
    unique_count = data[column_name].nunique()
    
    return unique_count

if __name__ == '__main__':
    # Check if the required arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <file_path> <column_index>")
        sys.exit(1)
    
    # Extract the arguments
    file_path = sys.argv[1]
    column_index = int(sys.argv[2])
    
    # Call the function to count unique elements
    unique_count = count_unique_elements(file_path, column_index)
    
    print(f"The number of unique elements in column at index {column_index} is: {unique_count}")
