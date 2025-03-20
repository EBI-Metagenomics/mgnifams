#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import os
import argparse
import csv
import pandas as pd

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--domtbl",
        required=True,
        metavar="FILE",
        type=str,
        help="TSV hmmsearch domtbl out results for filtering.",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the input folder with metadata files to get family-rep mapping and sizes.",
    )
    parser.add_argument(
        "-l",
        "--length_threshold",
        required=True,
        metavar="FLOAT",
        type=float,
        help="Minimum length percentage threshold of annotated domain (env) against query to keep.",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output txt file with the non-redundant family ids.",
    )
    return parser.parse_args(args)

def create_size_dict(folder_path):
    family_to_size = {}
    
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            base_name = os.path.splitext(filename)[0]  # Remove extension
            
            with open(file_path, newline='') as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    if len(row) < 6:
                        continue  # Skip malformed lines
                    
                    id_col, size, _, _, _, _ = row  # Extract columns
                    key = f"{base_name}_{id_col}"
                    family_to_size[key] = size
    
    return family_to_size

def create_rep_to_family_dict(folder_path, family_to_size, redundant_fam_names):
    representative_to_family = {}
    
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            base_name = os.path.splitext(filename)[0]  # Remove extension
            
            with open(file_path, newline='') as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    if len(row) < 6:
                        continue  # Skip malformed lines
                    
                    id_col, size, rep, region, length, sequence = row  # Extract columns
                    key = f"{rep}/{region}" if region != "'-" else rep
                    value = f"{base_name}_{id_col}"
                    
                    if key in representative_to_family:
                        existing_value = representative_to_family[key]
                        if int(family_to_size[value]) > int(family_to_size[existing_value]):
                            redundant_fam_names.add(existing_value)
                            representative_to_family[key] = value
                        else:
                            redundant_fam_names.add(value)
                    else:
                        representative_to_family[key] = value
    
    return representative_to_family

def remove_self_hits(domtbl_df):
    return domtbl_df[domtbl_df["target name"] != domtbl_df["query name"]]

def filter_by_length(domtbl_df, length_threshold):
    domtbl_df = domtbl_df[
        (domtbl_df["env to"] - domtbl_df["env from"] + 1) / domtbl_df["qlen"]
        >= length_threshold
    ]

    return domtbl_df

def identify_redundant_fams(domtbl, metadata, length_threshold, out_file):
    redundant_fam_names = set()

    family_to_size = create_size_dict(metadata)
    representative_to_family = create_rep_to_family_dict(metadata, family_to_size, redundant_fam_names)

    domtbl_df = pd.read_csv(
        domtbl, sep=r"\s+", comment="#", header=None, usecols=[0, 3, 5, 19, 20]
    ).rename(
        columns={
            0: "target name",
            3: "query name",
            5: "qlen",
            19: "env from",
            20: "env to",
        }
    )

    # Removing already found redundant families (those that had the exact same rep/region)
    domtbl_df = domtbl_df[~domtbl_df["query name"].isin(redundant_fam_names)]
    domtbl_df["target name"] = domtbl_df["target name"].map(representative_to_family)
    domtbl_df = remove_self_hits(domtbl_df)
    domtbl_df = filter_by_length(domtbl_df, length_threshold)
    domtbl_df = domtbl_df.drop(columns=["qlen", "env from", "env to"])
    domtbl_df["query size"] = domtbl_df["query name"].map(family_to_size)
    domtbl_df["target size"] = domtbl_df["target name"].map(family_to_size)

    # This method might deem redundant longer (super)families and keep shorter families with more members
    for _, row in domtbl_df.iterrows():
        if int(row["query size"]) < int(row["target size"]):
            redundant_fam_names.add(row["query name"])
        else:
            redundant_fam_names.add(row["target name"])
    
    with open(out_file, "w") as f:
        for item in sorted(redundant_fam_names):  # Sorting for consistent order
            f.write(f"{item}\n")

def main(args=None):
    args = parse_args(args)

    identify_redundant_fams(
        args.domtbl,
        args.metadata,
        args.length_threshold,
        args.out_file,
    )

if __name__ == "__main__":
    main()
