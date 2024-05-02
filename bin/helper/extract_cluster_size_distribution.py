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

def plot_hist(grouped_size_counts_df):
    # Separate values greater than 1000 and less than or equal to 1000
    above_1000 = grouped_size_counts_df[grouped_size_counts_df['Size'] > 1000]
    below_1000 = grouped_size_counts_df[grouped_size_counts_df['Size'] <= 1000]

    # Group values greater than 1000 into a single bin
    sum_above_1000 = above_1000['Count'].sum()

    # Calculate sums for every 50 units
    sums = []
    labels = []

    for i in range(0, 1000, 50):
        label_start = i + 1
        label_end = min(i + 50, 1000)
        label = f"{label_start}-{label_end}"

        if i == 950:
            sums.append(min(below_1000[below_1000['Size'] == i]['Count'].sum() + sum_above_1000, 1000))
        else:
            sums.append(min(below_1000[(below_1000['Size'] > i) & (below_1000['Size'] <= i + 50)]['Count'].sum(), 1000))

        labels.append(label)

    # Add the '>1000' label if there are counts above 1000
    if sum_above_1000 > 0:
        labels.append('>=1001')
        sums.append(min(sum_above_1000, 1000))

    # Create a new DataFrame with the combined data
    combined_df = pd.DataFrame({'Size': labels, 'Count': sums})

    # Plot the histogram
    plt.bar(labels, combined_df['Count'], width=0.8, edgecolor='black')
    plt.xlabel('Size')
    plt.ylabel('Count')
    plt.title('Histogram for 50 AA clusters ( >1000 Size Grouped and Count cut)')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

# Calc and Save
# file_path = '/nfs/production/rdf/metagenomics/users/vangelis/mgnifams/output/linclust_clusterss.tsv' # linclust result file
# # file_path = '/home/vangelis/Desktop/Projects/mgnifams-site/data/mmseqs_families_test.tsv' # local testing
# first_column_counts_df = count_first_column_elements(file_path)
# grouped_size_counts_df = group_and_count_by_size(first_column_counts_df)
# grouped_size_counts_df.to_csv('tmp/distr/grouped_size_counts.csv', index=False)
# # log10 also
# grouped_size_counts_log_df = grouped_size_counts_df.copy()
# grouped_size_counts_log_df['Count'] = np.log10(grouped_size_counts_log_df['Count'])
# grouped_size_counts_log_df.to_csv('tmp/distr/grouped_size_counts_log.csv', index=False)

# OR Load
grouped_size_counts_df = pd.read_csv('tmp/distr/grouped_size_counts_50_AA.csv')

# AND Plot
# plot_size_counts(grouped_size_counts_df)
plot_hist(grouped_size_counts_df)
