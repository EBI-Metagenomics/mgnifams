#!/usr/bin/env python3

import sys
import os
import polars as pl
import math

def parse_args():
    global arg_clusters, arg_checked_clusters, \
        arg_minimum_members, arg_num_cluster_chunks
        
    if not (len(sys.argv) == 5):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_clusters           = sys.argv[1]
    arg_checked_clusters   = sys.argv[2]
    arg_minimum_members    = int(sys.argv[3])
    arg_num_cluster_chunks = int(sys.argv[4])

def remove_clusters(df, removal_file):
    """
    Remove clusters from the DataFrame based on the cluster representatives listed in the removal file.
    
    Parameters:
    df (pl.DataFrame): The original Polars DataFrame.
    removal_file (str): Path to the TSV file containing cluster representative names to be removed.
    
    Returns:
    pl.DataFrame: The filtered DataFrame with specified clusters removed.
    """
    if (removal_file == '0'): # default code when no file is provided
        return df
    if os.path.getsize(removal_file) == 0: # if file provided, but empty
        return df

    # Read the removal file into a Polars DataFrame
    removal_df = pl.read_csv(removal_file, separator='\t', has_header=False, new_columns=["rep"])
    # Get a list of cluster representatives to be removed
    removal_list = removal_df["rep"].to_list()
    # Filter the original DataFrame to exclude the specified cluster representatives
    filtered_df = df.filter(~pl.col("rep").is_in(removal_list))

    return filtered_df

def filter_clusters_by_member_count(df, minimum_members):
    """
    Group by cluster representative and filter out clusters with fewer than the specified minimum members.

    Parameters:
    df (pl.DataFrame): The original Polars DataFrame.
    minimum_members (int): The minimum number of members required for a cluster to be retained.

    Returns:
    pl.DataFrame: The filtered DataFrame with clusters having fewer than the specified minimum members removed.
    """
    filtered_df = (
        df.group_by("rep")
        .agg(pl.col("mem").count().alias("mem_count"))
        .filter(pl.col("mem_count") >= arg_minimum_members)
        .join(df, on="rep", how="inner")
        .drop("mem_count")
    )

    return filtered_df

def split_and_write_chunks(df, num_chunks, out_dir):
    """
    Split the DataFrame into a specified number of chunks based on cluster representatives
    and write each chunk to a TSV file.

    Parameters:
    df (pl.DataFrame): The DataFrame to be split.
    num_chunks (int): The number of chunks to split the DataFrame into.
    """
    # Get unique cluster representatives
    unique_reps = df.select("rep").unique()
    # Determine the size of each chunk
    chunk_size = math.ceil(len(unique_reps) / num_chunks)
    # Split the unique cluster representatives into chunks
    for i in range(num_chunks):
        start_idx = i * chunk_size
        if i == num_chunks - 1:  # Last chunk
            end_idx = len(unique_reps)
        else:
            end_idx = start_idx + chunk_size
        chunk_reps = unique_reps[start_idx:end_idx].to_series().to_list()
        # Filter the original DataFrame for the current chunk of cluster representatives
        chunk_df = df.filter(pl.col("rep").is_in(chunk_reps))
        
        # Write the chunk to a TSV file
        if not chunk_df.is_empty():
            chunk_df.write_csv(os.path.join(out_dir, f"linclust_clusters_{i + 1}.tsv"), separator='\t', include_header=False)

def main():
    parse_args()

    out_dir = 'cluster_chunks'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    df = pl.read_csv(arg_clusters, separator='\t', has_header=False)
    df = df.rename({"column_1": "rep", "column_2": "mem"})
    
    filtered_df = remove_clusters(df, arg_checked_clusters)
    filtered_df = filter_clusters_by_member_count(filtered_df, arg_minimum_members)
    split_and_write_chunks(filtered_df, arg_num_cluster_chunks, out_dir)

if __name__ == "__main__":
    main()
