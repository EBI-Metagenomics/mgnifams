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
        "-d", "--domtbl", required=True, metavar="FILE", type=str,
        help="TSV hmmsearch domtbl out results for filtering.",
    )
    parser.add_argument(
        "-m", "--metadata", required=True, metavar="FOLDER", type=str,
        help="Name of the input folder with metadata files to get family-rep mapping and sizes.",
    )
    parser.add_argument(
        "-e", "--edgelist", required=True, metavar="FILE", type=str,
        help="Generated family TSV (rep\\tmember)).",
    )
    parser.add_argument(
        "-l", "--length_threshold", required=True, metavar="FLOAT", type=float,
        help="Minimum length percentage threshold of annotated domain (env) against query to keep.",
    )
    parser.add_argument(
        "-r", "--redundant_score_threshold", required=True, metavar="FLOAT", type=float,
        help="Minimum Jaccard Score between family proteins, to decide if one (the one with fewer members) will be considered redundant.",
    )
    parser.add_argument(
        "-s", "--similarity_score_threshold", required=True, metavar="FLOAT", type=float,
        help="Minimum Jaccard Score between family proteins, to decide if they will be annotated as similar.",
    )
    parser.add_argument(
        "-o", "--out_file", required=True, metavar="FILE", type=str,
        help="Name of the output txt file with the non-redundant family ids.",
    )
    parser.add_argument(
        "-c", "--similarity_csv", required=True, metavar="FILE", type=str,
        help="Name of the output csv file with the similar family ids.",
    )
    return parser.parse_args(args)

def create_size_dict(folder_path):
    family_to_size = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            base_name = os.path.splitext(filename)[0]
            with open(file_path, newline='') as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    if len(row) < 7:
                        continue
                    id_col, size, *_ = row
                    key = f"{base_name}_{id_col}"
                    family_to_size[key] = size
    return family_to_size

def create_rep_to_family_dict(folder_path, family_to_size, redundant_fam_names):
    representative_to_family = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            base_name = os.path.splitext(filename)[0]
            with open(file_path, newline='') as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    if len(row) < 6:
                        continue
                    id_col, _, rep, region, *_ = row
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
    return domtbl_df[
        (domtbl_df["env to"] - domtbl_df["env from"] + 1) / domtbl_df["qlen"]
        >= length_threshold
    ]

def jaccard_index(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union > 0 else 0.0

def write_similarity_csv(similarity_csv, similarity_rows):
    with open(similarity_csv, "w", newline='') as out_csv:
        # Write the comment lines manually
        out_csv.write(
            "# id: \"family_similarities\"\n"
            "# section_name: \"Inter-family Jaccard scores\"\n"
            "# description: \"Table of similar, but non-redundant protein families.\"\n"
            "# format: \"csv\"\n"
            "# plot_type: \"table\"\n"
            "Row,Family Id 1,Family Id 2,Jaccard Score\n"
        )

        if similarity_rows:
            for i, (fam1, fam2, score) in enumerate(similarity_rows, start=1):
                out_csv.write(f'{i},"{fam1}","{fam2}",{score}\n')
        else:
            pass

def identify_redundant_fams(domtbl, metadata, edgelist, length_threshold,
                            redundant_score_threshold, similarity_score_threshold,
                            out_file, similarity_csv):

    redundant_fam_names = set()
    family_to_size = create_size_dict(metadata)
    print("✅ Sizes dictionary created")
    representative_to_family = create_rep_to_family_dict(metadata, family_to_size, redundant_fam_names)
    print("✅ Rep to family dictionary created")

    domtbl_df = pd.read_csv(
        domtbl, sep=r"\s+", comment="#", header=None, usecols=[0, 3, 5, 19, 20]
    ).rename(columns={
        0: "target name",
        3: "query name",
        5: "qlen",
        19: "env from",
        20: "env to",
    })

    domtbl_df = domtbl_df[~domtbl_df["query name"].isin(redundant_fam_names)]
    domtbl_df["target name"] = domtbl_df["target name"].map(representative_to_family)
    domtbl_df = remove_self_hits(domtbl_df)
    domtbl_df = filter_by_length(domtbl_df, length_threshold)
    domtbl_df = domtbl_df.drop(columns=["qlen", "env from", "env to"])
    domtbl_df["query size"] = domtbl_df["query name"].map(family_to_size)
    domtbl_df["target size"] = domtbl_df["target name"].map(family_to_size)

    edge_df = pd.read_csv(edgelist, sep="\t", header=None, names=["family", "member"])
    print("✅ Families TSV pandas dataframe created")
    fam_to_proteins = edge_df.groupby("family")["member"].apply(list).to_dict()
    fam_to_base_proteins = {
        fam: set(m.split("/")[0] for m in members)
        for fam, members in fam_to_proteins.items()
    }
    print("✅ Family protein groups created and sliced")

    similarity_rows = []
    for _, row in domtbl_df.iterrows():
        fam1 = row["query name"]
        fam2 = row["target name"]

        set1 = fam_to_base_proteins.get(fam1, set())
        set2 = fam_to_base_proteins.get(fam2, set())
        score = jaccard_index(set1, set2)
        if score >= redundant_score_threshold:
            if int(row["query size"]) < int(row["target size"]):
                redundant_fam_names.add(fam1)
            elif int(row["query size"]) > int(row["target size"]):
                redundant_fam_names.add(fam2)
            else: # sizes equal, keep alphabetically first as non-redundant to avoid triangular deletions
                redundant_fam_names.add(max(fam1, fam2))
        elif score >= similarity_score_threshold:
            similarity_rows.append((fam1, fam2, score))
    print("✅ Jaccard scores calculated")

    write_similarity_csv(similarity_csv, similarity_rows)
    print("✅ Similarity CSV written")

    with open(out_file, "w") as f:
        for item in sorted(redundant_fam_names):
            f.write(f"{item}\n")
    print("✅ Redundant ids TXT written")

def main(args=None):
    args = parse_args(args)
    identify_redundant_fams(
        args.domtbl,
        args.metadata,
        args.edgelist,
        args.length_threshold,
        args.redundant_score_threshold,
        args.similarity_score_threshold,
        args.out_file,
        args.similarity_csv
    )

if __name__ == "__main__":
    main()
