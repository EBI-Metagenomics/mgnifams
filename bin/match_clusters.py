import sys
import pandas as pd

def filter_blastp_results(blastp_results_df):
    blastp_results_df["evalue"] = pd.to_numeric(blastp_results_df["evalue"], errors='coerce') # Convert 'evalue' column to numeric
    filtered_df = blastp_results_df[blastp_results_df["evalue"] <= 1e-7] # Filter rows with evalue <= 1e-7
    return filtered_df

def match_clusters(filtered_blastp_df, clusters_df):
    filtered_blastp_df['sseqid'] = filtered_blastp_df['sseqid'].astype(str)  # Convert 'sseqid' to string
    matched_rows = pd.merge(filtered_blastp_df, clusters_df, left_on='sseqid', right_on='Member', how='inner')
    matched_rows = matched_rows[["qseqid", "sseqid", "Representative"]]
    matched_rows.columns = ["anti_defence_protein_id", "matched_mgnify_id", "matched_cluster_rep_id"]
    return matched_rows

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py blastp_results clusters")
        sys.exit(1)

    blastp_results_file = sys.argv[1]
    clusters_file = sys.argv[2]

    # Read blastp_results file into a dataframe
    blastp_results_df = pd.read_csv(blastp_results_file, sep='\t', header=None)
    blastp_results_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    # Filter blastp results based on evalue
    filtered_blastp_df = filter_blastp_results(blastp_results_df)

    # Read clusters file into a dataframe
    clusters_df = pd.read_csv(clusters_file, sep='\t', header=None)
    clusters_df.columns = ["Representative", "Member"]

    # Match remaining sseqids with clusters
    matched_rows = match_clusters(filtered_blastp_df, clusters_df)
    matched_rows.to_csv("matched_cluster_reps.csv", index=False)
