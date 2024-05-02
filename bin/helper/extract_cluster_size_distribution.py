import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

def count_first_column_elements(file_path):
    print("count_first_column_elements\n")
    start_time = time.time()

    # Read the TSV file
    data = pd.read_csv(file_path, sep='\t', header=None)

    # Extract the first column
    first_column = data.iloc[:, 0]

    # Count the occurrences of each element
    counts = first_column.value_counts().reset_index()

    # Rename the columns
    counts.columns = ['Element', 'Size']

    print(str(time.time() - start_time) + "\n")

    return counts

def group_and_count_by_size(input_df):
    print("group_and_count_by_size\n")
    start_time = time.time()

    # Group by count and count occurrences
    size_counts = input_df.groupby('Size').size().reset_index(name='Count')
    size_counts['Count'] = size_counts['Count'].astype(int)
    size_counts = size_counts.sort_values(by='Count', ascending=False)

    print(str(time.time() - start_time) + "\n")

    return size_counts
    
def plot_size_counts(size_counts_df):
    print("plot_size_counts\n")
    start_time = time.time()

    plt.figure(figsize=(10, 6))
    sns.stripplot(x='Size', y='Count', data=size_counts_df, jitter=True)

    # Rotate x-axis labels for better readability if needed
    plt.xticks(rotation=45)

    # Add labels and title
    plt.xlabel('Size')
    plt.ylabel('Count')
    plt.title('Dot Plot of Size vs Count')

    print(str(time.time() - start_time) + "\n")

    # Show plot
    plt.show()

# Calc and Save
file_path = '/nfs/production/rdf/metagenomics/users/vangelis/mgnifams/output/linclust_clusterss.tsv' # linclust result file
# file_path = '/home/vangelis/Desktop/Projects/mgnifams-site/data/mmseqs_families_test.tsv' # local testing
first_column_counts_df = count_first_column_elements(file_path)
grouped_size_counts_df = group_and_count_by_size(first_column_counts_df)
grouped_size_counts_df.to_csv('tmp/distr/grouped_size_counts.csv', index=False)
# log10 also
grouped_size_counts_log_df = grouped_size_counts_df.copy()
grouped_size_counts_log_df['Count'] = np.log10(grouped_size_counts_log_df['Count'])
grouped_size_counts_log_df.to_csv('tmp/distr/grouped_size_counts_log.csv', index=False)

# OR Load
# grouped_size_counts_df = pd.read_csv('data/grouped_size_counts.csv')

# AND Plot
plot_size_counts(grouped_size_counts_df)
