#!/usr/bin/env python3
"""Parquet twin of query_mgnprotein_db.py for biomes: writes query_results/<family>.tsv with one line per distinct
member protein, `mgyp<TAB><cluster_rep_biomes>` (`;`-joined biome lineages, empty if none; parse_biomes.py reads
column 1 without a mapping), from the family TSVs and the MGnify clusters parquet. Members must be cluster
representatives: their own biomes are their row's cluster_rep_biomes, as the MGnify proteins DB `b` was in 1.0.
"""

import argparse
import os
from collections import defaultdict

import pyarrow.dataset as ds
from parse_domains import extract_mgyp


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", required=True, nargs="+", help="family<TAB>member TSV(s)")
    parser.add_argument("--clusters", required=True, help="MGnify clusters parquet (cluster_rep, cluster_rep_biomes)")
    parser.add_argument("--output_dir", default="query_results")
    a = parser.parse_args(args)

    families = defaultdict(set)
    for tsv in a.families:
        with open(tsv) as fh:
            for line in fh:
                family, member = line.rstrip("\n").split("\t")
                families[family].add(int(extract_mgyp(member)))  # non-numeric protein ids fail here

    ids = sorted({p for members in families.values() for p in members})
    table = ds.dataset(a.clusters, format="parquet").to_table(
        columns=["cluster_rep", "cluster_rep_biomes"], filter=ds.field("cluster_rep").isin(ids)
    )
    biomes = dict(zip(table.column("cluster_rep").to_pylist(), table.column("cluster_rep_biomes").to_pylist()))
    missing = [p for p in ids if p not in biomes]
    if missing:
        raise SystemExit(f"{len(missing)} members are not cluster_rep rows of {a.clusters}, e.g. {missing[:10]}")

    os.makedirs(a.output_dir, exist_ok=True)
    for family, members in families.items():
        with open(os.path.join(a.output_dir, f"{family}.tsv"), "w") as out:
            for protein in sorted(members):
                out.write(f"{protein}\t{biomes[protein] or ''}\n")


if __name__ == "__main__":
    main()
