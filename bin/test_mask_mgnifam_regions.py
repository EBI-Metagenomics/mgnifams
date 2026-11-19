#!/usr/bin/env python3
"""Self-check of mask_mgnifam_regions.py. Run: python3 bin/test_mask_mgnifam_regions.py"""

import random
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN))
from mask_mgnifam_regions import fasta_id, mask_record, region  # noqa: E402
from slice_protein import sliceProtein  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def pfam_pass(protein, seq, pfam, min_len):
    """What the Pfam slicer writes for one protein (extract_unannotated_parquet_slices.py)."""
    if pfam:
        return sliceProtein(protein, seq, {"p": [["PF", s, e] for s, e in pfam]}, min_len)
    return [f">{protein}\n{seq}\n"] if len(seq) >= min_len else []


def one_pass(protein, seq, pfam, mgnifam, min_len):
    """Pfam and MGnifam regions removed together; a whole-protein MGnifam row (end 0) covers the protein."""
    regions = [["p", s, e] for s, e in pfam] + [["m", s, e or len(seq)] for s, e in mgnifam]
    if regions:
        return sliceProtein(protein, seq, {"p": regions}, min_len)
    return [f">{protein}\n{seq}\n"] if len(seq) >= min_len else []


def two_pass(protein, seq, pfam, mgnifam, min_len):
    out = []
    for record in pfam_pass(protein, seq, pfam, min_len):
        header, sequence = record[1:].split("\n")[:2]
        out += mask_record(header, sequence, mgnifam, min_len)
    return out


def random_regions(rng, length, n):
    regs = []
    for _ in range(n):
        s = rng.randint(1, length)
        regs.append((s, rng.randint(s, length)))
    return regs


# (a) one pass = two passes, on random proteins
rng = random.Random(20261119)
for case in range(20000):
    length = rng.randint(1, 600)
    seq = "".join(rng.choice(AA) for _ in range(length))
    pfam = random_regions(rng, length, rng.randint(0, 3))
    mgnifam = random_regions(rng, length, rng.randint(0, 4))
    if rng.random() < 0.05:
        mgnifam.append((1, 0))
    min_len = rng.choice((1, 50, 100))
    want = sorted(one_pass(7, seq, pfam, mgnifam, min_len))
    got = sorted(two_pass(7, seq, pfam, mgnifam, min_len))
    assert got == want, (case, length, pfam, mgnifam, min_len, got, want)

seq = AA * 5  # 100 aa
# (b) a region partly outside a slice is clipped; one wholly outside leaves the record unchanged
assert mask_record("9_21_60", seq[20:60], [(1, 30)], 5) == [f">9_31_60\n{seq[30:60]}\n"]
assert mask_record("9_21_60", seq[20:60], [(55, 90)], 5) == [f">9_21_54\n{seq[20:54]}\n"]
assert mask_record("9_21_60", seq[20:60], [(61, 90), (1, 20)], 5) == [f">9_21_60\n{seq[20:60]}\n"]
# (c) a whole-protein row (end 0) masks whole and sliced records
assert mask_record("9", seq, [(1, 0)], 1) == []
assert mask_record("9_21_60", seq[20:60], [(1, 0)], 1) == []
# (d) overlapping and repeated regions
assert mask_record("9", seq, [(10, 30), (25, 40), (25, 40)], 5) == [f">9_1_9\n{seq[:9]}\n", f">9_41_100\n{seq[40:]}\n"]
# a piece shorter than the minimum is dropped
assert mask_record("9", seq, [(4, 100)], 5) == []
# (e) padded 1.0 names, bare proteins
assert region("12\t10196381955/5-330   \n") == (10196381955, 5, 330)
assert region("12\t10196381955   \n") == (10196381955, 1, 0)
# (f) a record with no region keeps its header and sequence
assert mask_record("9_21_60", seq[20:60], [], 5) == [f">9_21_60\n{seq[20:60]}\n"]
# (g) both spellings give the same id, so the same bucket
assert fasta_id(">10196381955_5_330") == fasta_id(">10196381955") == region("1\t10196381955/5-330")[0] == 10196381955
# (h) non-integer ids and broken headers stop the run
for bad in (
    lambda: fasta_id(">MGYP000123"),
    lambda: region("1\tMGYP000123/1-5"),
    lambda: mask_record("9_5", "AAAA", [], 1),
):
    try:
        bad()
    except SystemExit:
        pass
    else:
        raise AssertionError("expected SystemExit")
try:
    mask_record("9_1_10", "AAAA", [], 1)  # header span differs from the sequence
except SystemExit:
    pass
else:
    raise AssertionError("expected SystemExit")
# (i) sliceProtein's default start is 1
assert sliceProtein(9, seq, {"m": [["m", 10, 30]]}, 5) == sliceProtein(9, seq, {"m": [["m", 10, 30]]}, 5, start=1)

# CLI: split into buckets, mask each bucket, joined output = masking without buckets
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    proteins = {p: "".join(rng.choice(AA) for _ in range(rng.randint(120, 400))) for p in range(100, 140)}
    with open(d / "in.fa", "w") as fh:
        for p, s in proteins.items():
            if p % 3:  # sliced record, wrapped over two lines
                fh.write(f">{p}_11_{len(s)}\n{s[10:60]}\n{s[60:]}\n")
            else:
                fh.write(f">{p}\n{s}\n")
    rows = {p: random_regions(rng, len(s), p % 4) for p, s in proteins.items()}
    (d / "a.tsv").write_text("".join(f"f1\t{p}/{s}-{e}\n" for p, regs in rows.items() for s, e in regs if p % 2))
    (d / "b.tsv").write_text(
        "".join(f"f2\t{p}/{s}-{e}  \n" for p, regs in rows.items() for s, e in regs if not p % 2) + "f3\t105\n"
    )
    rows[105] = rows[105] + [(1, 0)]
    run = [sys.executable, BIN / "mask_mgnifam_regions.py"]
    subprocess.run([*run, "split", "--kind", "fasta", "--buckets", "7", "--out", d / "b", d / "in.fa"], check=True)
    subprocess.run(
        [*run, "split", "--kind", "regions", "--buckets", "7", "--out", d / "b", d / "a.tsv", d / "b.tsv"], check=True
    )
    assert sorted(f.name for f in (d / "b").iterdir()) == sorted(
        [f"{i:03d}.fa" for i in range(7)] + [f"{i:03d}.tsv" for i in range(7)]
    )
    got = []
    for i in range(7):
        out = d / f"{i:03d}.masked.fa"
        subprocess.run(
            [
                *run,
                "mask",
                "--fasta",
                d / "b" / f"{i:03d}.fa",
                "--regions",
                d / "b" / f"{i:03d}.tsv",
                "--min_sequence_length",
                "30",
                "--out",
                out,
            ],
            check=True,
        )
        got += [">" + r for r in out.read_text().split(">") if r]
    want = []
    for p, s in proteins.items():
        start = 11 if p % 3 else 1
        header = f"{p}_11_{len(s)}" if p % 3 else str(p)
        want += mask_record(header, s[start - 1 :], rows[p], 30)
    assert sorted(got) == sorted(want), (len(got), len(want))
    assert not any(r.startswith(">105") for r in got)
    bad = subprocess.run(
        [*run, "split", "--kind", "fasta", "--buckets", "7", "--out", d / "x", d / "a.tsv"],
        capture_output=True,
        text=True,
    )
    assert bad.returncode != 0
print("test_mask_mgnifam_regions: OK")
