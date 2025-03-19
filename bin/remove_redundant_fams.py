#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import os
import sys
import pandas as pd
import argparse
import shutil


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
        "-f",
        "--fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the input file file with the pre-filtered reps fasta.",
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


def remove_self_hits(domtbl_df, representative_to_family):
    domtbl_df["target name"] = domtbl_df["target name"].map(representative_to_family)
    domtbl_df = domtbl_df[domtbl_df["target name"] != domtbl_df["query name"]]

    return domtbl_df


def filter_by_length(domtbl_df, length_threshold):
    domtbl_df = domtbl_df[
        (domtbl_df["env to"] - domtbl_df["env from"] + 1) / domtbl_df["qlen"]
        >= length_threshold
    ]

    return domtbl_df


def remove_redundant_fams(domtbl, fasta, metadata, length_threshold, out_folder):
    pass
    # TODO create dictionary of reps.fasta of seq name -> family id
    # mapping_df = pd.read_csv(
    #     mapping, comment="#", usecols=["Family Id", "Size", "Representative Id"]
    # )
    # domtbl_df = pd.read_csv(
    #     domtbl, sep=r"\s+", comment="#", header=None, usecols=[0, 3, 5, 19, 20]
    # ).rename(
    #     columns={
    #         0: "target name",
    #         3: "query name",
    #         5: "qlen",
    #         19: "env from",
    #         20: "env to",
    #     }
    # )

    # representative_to_family = dict(
    #     zip(mapping_df["Representative Id"], mapping_df["Family Id"])
    # )
    # family_to_size = dict(zip(mapping_df["Family Id"], mapping_df["Size"]))

    # domtbl_df = remove_self_hits(domtbl_df, representative_to_family)
    # domtbl_df = filter_by_length(domtbl_df, length_threshold)
    # domtbl_df = domtbl_df.drop(columns=["qlen", "env from", "env to"])
    # domtbl_df["query size"] = domtbl_df["query name"].map(family_to_size)
    # domtbl_df["target size"] = domtbl_df["target name"].map(family_to_size)

    # redundant_fam_names = set()
    # for _, row in domtbl_df.iterrows():
    #     if row["query size"] < row["target size"]:
    #         redundant_fam_names.add(row["query name"])
    #     else:
    #         redundant_fam_names.add(row["target name"])

    # for file_name in os.listdir(fasta_folder):
    #     base_name = os.path.basename(file_name).split(".")[0]
    #     if base_name not in redundant_fam_names:
    #         source_file = os.path.join(fasta_folder, file_name)
    #         destination_file = os.path.join(out_folder, file_name)

    #         # Check if it is a file (not a directory) and copy
    #         if os.path.isfile(source_file):
    #             shutil.copy2(source_file, destination_file)


def main(args=None):
    args = parse_args(args)

    remove_redundant_fams(
        args.domtbl,
        args.fasta,
        args.metadata,
        args.length_threshold,
        args.out_file,
    )


if __name__ == "__main__":
    sys.exit(main())
