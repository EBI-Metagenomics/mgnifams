#!/usr/bin/env python3
"""Parquet twin of query_mgnprotein_db.py for biomes: writes query_results/<family>.tsv with one line per distinct
member protein, `mgyp<TAB><cluster_rep_biomes>` (`;`-joined biome lineages, empty if none; parse_biomes.py reads
column 1 without a mapping), from the family TSVs and the MGnify clusters parquet. Members must be cluster
representatives: their own biomes are their row's cluster_rep_biomes, as the MGnify proteins DB `b` was in 1.0.
"""

import argparse
import os

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
from parse_domains import load_families


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", required=True, nargs="+", help="family<TAB>member TSV(s)")
    parser.add_argument("--clusters", required=True, help="MGnify clusters parquet (cluster_rep, cluster_rep_biomes)")
    parser.add_argument("--output_dir", default="query_results")
    a = parser.parse_args(args)

    families = load_families(a.families)
    # All members, family after family (sorted within each), as one int64 array: no Python object per member
    members = pa.chunked_array(
        [pa.array(sorted(m), pa.int64()) for m in families.values()], pa.int64()
    ).combine_chunks()

    # Read dictionary-encoded: the column holds few distinct lineage strings, one Python string each
    fmt = ds.ParquetFileFormat(read_options={"dictionary_columns": ["cluster_rep_biomes"]})
    table = (
        ds.dataset(a.clusters, format=fmt)
        .to_table(
            columns=["cluster_rep", "cluster_rep_biomes"], filter=ds.field("cluster_rep").isin(pc.unique(members))
        )
        .unify_dictionaries()
        .combine_chunks()
    )
    positions = pc.index_in(members, value_set=table.column("cluster_rep"))
    missing = pc.unique(pc.filter(members, pc.is_null(positions)))
    if len(missing):
        raise SystemExit(
            f"{len(missing)} members are not cluster_rep rows of {a.clusters}, e.g. {missing[:10].to_pylist()}"
        )

    biomes = table.column("cluster_rep_biomes").combine_chunks()
    lineages = biomes.dictionary.to_pylist() + [""]  # last: members without biomes (null)
    codes = pc.take(biomes.indices, positions).fill_null(len(lineages) - 1)

    os.makedirs(a.output_dir, exist_ok=True)
    offset = 0
    for family, family_members in families.items():
        n = len(family_members)
        with open(os.path.join(a.output_dir, f"{family}.tsv"), "w") as out:
            for protein, code in zip(members.slice(offset, n).to_pylist(), codes.slice(offset, n).to_pylist()):
                out.write(f"{protein}\t{lineages[code]}\n")
        offset += n


if __name__ == "__main__":
    main()
