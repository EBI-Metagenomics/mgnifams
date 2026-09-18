#!/usr/bin/env python3
"""Equivalence check: parquet slicer == CSV slicer on the same proteins. Run: python3 bin/test_extract_unannotated_parquet_slices.py"""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parent
MIN = 10
seq = "ACDEFGHIKLMNPQRSTVWY" * 3  # 60 aa
# protein_id -> Pfam [env_from, env_to] regions
proteins = {
    1: [[5, 20], [18, 30]],  # overlapping
    2: [[11, 20], [21, 30]],  # adjacent
    3: [[1, 45]],  # slice runs to the protein end
    4: [],  # no Pfam rows -> whole protein
    5: [[20, 40]],  # both flanks >= MIN
    6: [[2, 55]],  # nothing left
    7: [],  # shorter than MIN -> dropped
}
seqs = {p: (seq[:8] if p == 7 else seq) for p in proteins}


def records(fasta):
    text = Path(fasta).read_text()
    return sorted(">" + r for r in text.split(">") if r)


with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    with open(d / "in.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mgyp", "sequence", "full_length", "cluster_size", "metadata"])
        for p, regs in proteins.items():
            meta = {"p": [["PF00001", 1e-5, 30.0, 1, 10, s, e] for s, e in regs]} if regs else {"b": [[1, 1]]}
            w.writerow([p, seqs[p], "true", 1, json.dumps(meta)])
    subprocess.run(
        [
            sys.executable,
            BIN / "extract_unannotated_slices.py",
            "-i",
            d / "in.csv",
            "-o",
            d / "csv.faa",
            "-l",
            str(MIN),
        ],
        check=True,
    )

    ids = list(proteins)
    pq.write_table(
        pa.table({"protein_id": pa.array(ids, pa.int64()), "sequence": [seqs[p] for p in ids]}),
        d / "seqs.parquet",
        row_group_size=2,  # 4 row groups
    )
    rows = [(p, s, e) for p, regs in proteins.items() for s, e in regs]
    pq.write_table(
        pa.table({k: pa.array(v, pa.int64()) for k, v in zip(["protein_id", "env_from", "env_to"], zip(*rows))}),
        d / "pfam.parquet",
    )
    parquet_records = []
    for n_chunks in (1, 3):
        parquet_records = []
        for i in range(n_chunks):
            out = d / f"chunk_{i}.faa"
            subprocess.run(
                [
                    sys.executable,
                    BIN / "extract_unannotated_parquet_slices.py",
                    "--sequences",
                    d / "seqs.parquet",
                    "--pfam",
                    d / "pfam.parquet",
                    "--chunk_index",
                    str(i),
                    "--n_chunks",
                    str(n_chunks),
                    "--min_sequence_length",
                    str(MIN),
                    "--output_file",
                    out,
                ],
                check=True,
            )
            parquet_records += records(out)
        assert sorted(parquet_records) == records(d / "csv.faa"), (n_chunks, parquet_records, records(d / "csv.faa"))

    # more chunks than row groups -> empty output, no error
    subprocess.run(
        [
            sys.executable,
            BIN / "extract_unannotated_parquet_slices.py",
            "--sequences",
            d / "seqs.parquet",
            "--pfam",
            d / "pfam.parquet",
            "--chunk_index",
            "5",
            "--n_chunks",
            "6",
            "--min_sequence_length",
            str(MIN),
            "--output_file",
            d / "empty.faa",
        ],
        check=True,
    )
    assert (d / "empty.faa").read_text() == ""
    assert any(r.startswith(">3_46_60") for r in parquet_records) and any(r.startswith(">4\n") for r in parquet_records)
print("ok")
