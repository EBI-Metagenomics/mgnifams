#!/usr/bin/env python3
"""Self-check for update_sqlite_blobs_staged.py. Run: python3 bin/test_update_sqlite_blobs_staged.py"""

import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = Path(__file__).resolve().parent
SCHEMA = BIN.parent / "assets/data/db_schema.sqlite"


def run(d, *extra):
    cmd = [sys.executable, BIN / "update_sqlite_blobs_staged.py", "--db", d / "db.sqlite3", "--ids", d / "ids.txt"]
    cmd += ["--cif_dir", d / "cif", "--s4pred_dir", d / "s4pred", "--domain_dir", d / "domain"]
    return subprocess.run([*cmd, "--required", "cif,s4pred,domain", *extra], capture_output=True, text=True)


with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    con = sqlite3.connect(d / "db.sqlite3")
    con.executescript(SCHEMA.read_text())
    for i in (1, 2):
        con.execute(
            "INSERT INTO mgnifam (id, full_size, protein_rep, rep_region, rep_length, rep_sequence, plddt, ptm,"
            " helix_percent, strand_percent, coil_percent) VALUES (?, 5, 7, '1-9', 9, 'ACDEFGHIK', 80, 0.5, 1, 2, 97)",
            (i,),
        )
    con.commit()
    con.close()
    (d / "ids.txt").write_text("1\n2\n")
    for sub, ext in (("cif", "cif"), ("s4pred/json", "json"), ("domain", "json")):
        (d / sub).mkdir(parents=True)
        for i in (1, 2):
            (d / sub / f"{i}.{ext}").write_text(f"{sub} {i}")

    (d / "cif/2.cif").unlink()
    r = run(d)
    assert r.returncode == 1 and "['2']" in r.stderr, r.stderr
    blobs = sqlite3.connect(d / "db.sqlite3").execute("SELECT count(cif_blob) FROM mgnifam").fetchone()[0]
    assert blobs == 0, "failed load must roll back"

    (d / "cif/2.cif").write_text("")
    assert "empty file" in run(d).stderr

    (d / "cif/2.cif").write_text("cif 2")
    (d / "ids.txt").write_text("1\n2\n3\n")
    assert "missing ['3']" in run(d).stderr
    (d / "ids.txt").write_text("1\n2\n")

    r = run(d, "--info", "tm_computed=false", "child_tables=pfams,funfams,folds")
    assert r.returncode == 0, r.stderr
    con = sqlite3.connect(d / "db.sqlite3")
    assert con.execute("SELECT count(*) FROM mgnifam WHERE length(s4pred_blob) > 0").fetchone()[0] == 2
    assert dict(con.execute("SELECT key, value FROM update_info"))["child_tables"] == "pfams,funfams,folds"
print("ok")
