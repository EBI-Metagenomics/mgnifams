#!/usr/bin/env python3
"""Parquet twin of extract_unannotated_slices.py: slice the regions of MGnify proteins not covered by known Pfam
domains (same sliceProtein logic), reading the sequence and Pfam parquet files joined on protein_id.

Chunk `chunk_index` of `n_chunks` processes the sequence row groups rg with rg % n_chunks == chunk_index.
"""

import argparse

import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from extract_unannotated_slices import sliceProtein


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True, help="Parquet with protein_id, sequence")
    parser.add_argument("--pfam", required=True, help="Parquet with protein_id, env_from, env_to (one row per domain)")
    parser.add_argument("--chunk_index", required=True, type=int)
    parser.add_argument("--n_chunks", required=True, type=int)
    parser.add_argument("--min_sequence_length", required=True, type=int)
    parser.add_argument("--output_file", required=True, help="Output FASTA")
    return parser.parse_args(args)


def main(args=None):
    a = parse_args(args)
    if a.n_chunks < 1 or not 0 <= a.chunk_index < a.n_chunks:
        raise SystemExit(f"invalid chunk {a.chunk_index} of {a.n_chunks}")

    sequences = pq.ParquetFile(a.sequences)
    pfam = ds.dataset(a.pfam, format="parquet")
    pid = ds.field("protein_id")

    with open(a.output_file, "w") as out:
        for rg in range(a.chunk_index, sequences.num_row_groups, a.n_chunks):
            table = sequences.read_row_group(rg, columns=["protein_id", "sequence"])
            if table.num_rows == 0:
                continue
            ids = table.column("protein_id")
            bounds = pc.min_max(ids)
            # ponytail: min/max prunes Pfam row groups only if that file is sorted by protein_id
            hits = pfam.to_table(
                columns=["protein_id", "env_from", "env_to"],
                filter=(pid >= bounds["min"]) & (pid <= bounds["max"]) & pid.isin(ids),
            )
            regions = {}
            for protein, env_from, env_to in zip(*(hits.column(c).to_pylist() for c in hits.column_names)):
                regions.setdefault(protein, []).append([env_from, env_to])

            for protein, sequence in zip(ids.to_pylist(), table.column("sequence").to_pylist()):
                if protein in regions:
                    out.writelines(sliceProtein(protein, sequence, {"p": regions[protein]}, a.min_sequence_length))
                elif len(sequence) >= a.min_sequence_length:
                    out.write(f">{protein}\n{sequence}\n")


if __name__ == "__main__":
    main()
