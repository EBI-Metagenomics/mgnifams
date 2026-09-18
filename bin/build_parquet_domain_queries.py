#!/usr/bin/env python3
"""Parquet twin of query_mgnprotein_db.py for domain architectures: writes query_results/<family>.tsv with one
line per distinct member protein, `mgyp<TAB><TAB>[['PF%05d', i_evalue, score, hmm_from, hmm_to, env_from, env_to], ...]`
(the biome column stays empty; parse_domains.py reads columns 0 and 2), from the refined families TSV and the MGnify
Pfam parquet file.
"""

import argparse
import os
from collections import defaultdict

import pyarrow.dataset as ds
from parse_domains import extract_mgyp

COLUMNS = ["protein_id", "pfam_accession", "i_evalue", "score", "hmm_from", "hmm_to", "env_from", "env_to"]


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refined_families", required=True, help="family<TAB>member TSV")
    parser.add_argument("--pfam", required=True, help="MGnify Pfam hits parquet")
    parser.add_argument("--output_dir", default="query_results")
    a = parser.parse_args(args)

    families = defaultdict(set)
    with open(a.refined_families) as fh:
        for line in fh:
            family, member = line.rstrip("\n").split("\t")
            families[family].add(int(extract_mgyp(member)))  # non-numeric protein ids fail here

    ids = sorted({p for members in families.values() for p in members})
    pfams = defaultdict(list)
    table = ds.dataset(a.pfam, format="parquet").to_table(columns=COLUMNS, filter=ds.field("protein_id").isin(ids))
    for protein, acc, i_evalue, score, *coords in zip(*(table.column(c).to_pylist() for c in COLUMNS)):
        pfams[protein].append(["PF%05d" % acc, i_evalue, score, *coords])

    os.makedirs(a.output_dir, exist_ok=True)
    for family, members in families.items():
        with open(os.path.join(a.output_dir, f"{family}.tsv"), "w") as out:
            for protein in sorted(members):
                out.write(f"{protein}\t\t{pfams[protein] if protein in pfams else ''}\n")


if __name__ == "__main__":
    main()
