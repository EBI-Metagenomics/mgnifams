#!/usr/bin/env python3
"""Redundant families from mgnifam compute_clans edges.

A pair is redundant when both its protein_jaccard (shared MGYP accessions) and its residue_jaccard reach their
thresholds. Families are visited from largest to smallest (full_msa_size; equal sizes: the name that sorts first
comes first); each family not yet dropped is kept and drops its redundant neighbours. On a chain A > B > C with
A~B and B~C this keeps A and C.

Values are compared as written by compute_clans (six decimals), so a pair within 5e-7 of a threshold may fall on
either side of it.
"""

import argparse
import csv
import os
from fractions import Fraction


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--edges", required=True, help="compute_clans <prefix>_clan_edges.csv")
    parser.add_argument("--metadata", required=True, nargs="+", help="mgnifam <chunk>_metadata.csv files")
    parser.add_argument("--score_threshold", required=True, help="minimum protein_jaccard")
    parser.add_argument("--residue_threshold", required=True, help="minimum residue_jaccard")
    parser.add_argument("--out", required=True, help="redundant family names, one per line")
    return parser.parse_args(args)


def read_sizes(paths):
    """`<chunk>_<family_id>` -> full_msa_size, the key pooling uses."""
    sizes = {}
    for path in paths:
        name = os.path.basename(path)
        if not name.endswith("_metadata.csv"):
            raise ValueError(f"{name}: expected a name ending in _metadata.csv")
        chunk = name.removesuffix("_metadata.csv")
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                sizes[f"{chunk}_{row['family_id']}"] = int(row["full_msa_size"])
    return sizes


def redundant_pairs(edges_path, score_threshold, residue_threshold):
    with open(edges_path, newline="") as handle:
        for row in csv.DictReader(handle):
            if (
                Fraction(row["protein_jaccard"]) >= score_threshold
                and Fraction(row["residue_jaccard"]) >= residue_threshold
            ):
                yield row["family_1"], row["family_2"]


def select_redundant(pairs, sizes):
    neighbours = {}
    for a, b in pairs:
        for family in (a, b):
            if family not in sizes:
                raise ValueError(f"{family}: in the edges but in no metadata file")
        neighbours.setdefault(a, set()).add(b)
        neighbours.setdefault(b, set()).add(a)
    kept, dropped = set(), set()
    for family in sorted(neighbours, key=lambda name: (-sizes[name], name)):
        if family in dropped:
            continue
        kept.add(family)
        dropped.update(neighbours[family] - kept)
    return sorted(dropped)


def main(args=None):
    args = parse_args(args)
    sizes = read_sizes(args.metadata)
    pairs = redundant_pairs(args.edges, Fraction(args.score_threshold), Fraction(args.residue_threshold))
    redundant = select_redundant(pairs, sizes)
    with open(args.out, "w") as handle:
        handle.writelines(f"{name}\n" for name in redundant)
    print(f"{len(redundant)} redundant families of {len(sizes)}")


if __name__ == "__main__":
    main()
