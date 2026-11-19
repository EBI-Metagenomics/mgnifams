#!/usr/bin/env python3
"""Check: parquet slicer == sliceProtein (bin/slice_protein.py) on the same proteins.
Run: python3 bin/test_extract_unannotated_parquet_slices.py"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from slice_protein import sliceProtein  # noqa: E402

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


# Output of the removed CSV slicer (extract_unannotated_slices.py, pipeline 3.x) on these proteins: sliceProtein
# with the default start must keep it byte for byte
GOLDEN = (
    ">1_31_60\nMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n>2_1_10\nACDEFGHIKL\n>2_31_60\nMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
    ">3_46_60\nGHIKLMNPQRSTVWY\n>4\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
    ">5_1_19\nACDEFGHIKLMNPQRSTVW\n>5_41_60\nACDEFGHIKLMNPQRSTVWY\n"
)
expected = []
for p, regs in proteins.items():
    if regs:
        expected += sliceProtein(p, seqs[p], {"p": [["PF00001", 1e-5, 30.0, 1, 10, s, e] for s, e in regs]}, MIN)
    elif len(seqs[p]) >= MIN:
        expected.append(f">{p}\n{seqs[p]}\n")
assert "".join(expected) == GOLDEN, expected
expected = sorted(expected)
# start = protein coordinate of sequence[0]: a slice keeps protein coordinates in its headers
assert sliceProtein(9, "A" * 30, {"m": [[11, 15]]}, 5, start=101) == [
    ">9_101_110\nAAAAAAAAAA\n",
    ">9_116_130\nAAAAAAAAAAAAAAA\n",
]

with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    # protein 8 is a cluster member, not a representative: the parquet slicer must drop it
    ids = [*proteins, 8]
    pq.write_table(
        pa.table({"protein_id": pa.array(ids, pa.int64()), "sequence": [seqs.get(p, seq) for p in ids]}),
        d / "seqs.parquet",
        row_group_size=2,  # 4 row groups
    )
    pq.write_table(pa.table({"cluster_rep": pa.array(list(proteins), pa.int64())}), d / "clusters.parquet")
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
                    "--clusters",
                    d / "clusters.parquet",
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
        assert sorted(parquet_records) == expected, (n_chunks, parquet_records, expected)

    # more chunks than row groups -> empty output, no error
    subprocess.run(
        [
            sys.executable,
            BIN / "extract_unannotated_parquet_slices.py",
            "--sequences",
            d / "seqs.parquet",
            "--clusters",
            d / "clusters.parquet",
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
