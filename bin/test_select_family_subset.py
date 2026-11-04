#!/usr/bin/env python3
"""Self-check for select_family_subset.py. Run: python3 bin/test_select_family_subset.py"""

import contextlib
import gzip
import io
import tempfile
from pathlib import Path

from select_family_subset import main

HEADER = (
    "family_id\tstatus\tnot_updated_reason\tfirst_release\tmodel_release\tmodel_proteins_release\t"
    "members_release\tmembers_proteins_release\tfull_size\tseed_size\trep_id\trep_region"
)


def setup(d, statuses, pdbs=("1", "2", "3", "9", "3_rep")):
    with gzip.open(d / "families.tsv.gz", "wt") as fh:
        fh.write(HEADER + "\n")
        for fid, status in statuses:
            fh.write(f"{fid}\t{status}\t\t1.0\t1.0\t2024_04\t1.0\t2024_04\t10\t5\tMGYP1\t1-10\n")
    (d / "reps.fasta").write_text(">1\nAAAA\n>2 desc text\nCCCC\nCC\n>3\nDDDD\n>9\nEEEE\n")
    (d / "structures").mkdir()
    for name in pdbs:
        (d / "structures" / f"{name}.pdb").write_text(f"MODEL {name}\n")


def run(d):
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        try:
            code = main(
                [
                    "--families", str(d / "families.tsv.gz"),
                    "--reps", str(d / "reps.fasta"),
                    "--structures", str(d / "structures"),
                    "--status", "not_updated",
                    "--out_fasta", str(d / "family_ids.fasta"),
                    "--out_structures", str(d / "pdb"),
                ]
            )
        except SystemExit as e:
            code = e.code
    return code or 0, err.getvalue()


def headers(d):
    return [line for line in (d / "family_ids.fasta").read_text().splitlines() if line.startswith(">")]


STATUSES = [("1", "updated"), ("2", "not_updated"), ("3", "not_updated"), ("4", "new")]

# selects_not_updated
with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    setup(d, STATUSES)
    code, err = run(d)
    assert code == 0, err
    assert headers(d) == [">2", ">3"], headers(d)
    assert sorted(p.name for p in (d / "pdb").iterdir()) == ["2.pdb", "3.pdb"]
    # header_first_token: description dropped, multi-line sequence kept
    assert (d / "family_ids.fasta").read_text() == ">2\nCCCC\nCC\n>3\nDDDD\n"

# missing_ids_fail
with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    setup(d, STATUSES, pdbs=("1", "2", "9", "3_rep"))
    code, err = run(d)
    assert code == 1, (code, err)
    assert "missing from structures: ['3']" in err, err

with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    setup(d, STATUSES + [("5", "not_updated")])
    code, err = run(d)
    assert code == 1 and "missing from reps: ['5']" in err, (code, err)

# zero_selected: exits 0 with empty outputs
with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    setup(d, [("1", "updated"), ("2", "updated")])
    code, err = run(d)
    assert code == 0, err
    assert (d / "family_ids.fasta").read_text() == ""
    assert list((d / "pdb").iterdir()) == []

print("ok")
