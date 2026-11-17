#!/usr/bin/env python3
"""Self-check for bin/load_release_manifest.py, with the repository's FTP manifests.
Run: python3 bin/test_load_release_manifest.py"""

import json
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = Path(__file__).resolve().parent
ROOT = BIN.parent
MANIFESTS = [str(ROOT / f"ftp/releases/{r}/release_manifest.json") for r in ("1.0", "1.1")]


def make_db(path, schema, releases):
    con = sqlite3.connect(path)
    con.executescript(schema)
    con.executemany(
        "INSERT INTO release_history (release, date, type, mgnify_proteins_release) VALUES (?, ?, 'update', '2026_07')",
        releases,
    )
    con.commit()
    con.close()


def run(db, *manifests):
    return subprocess.run(
        [sys.executable, str(BIN / "load_release_manifest.py"), "--db", db, "--manifests", *manifests],
        capture_output=True,
        text=True,
    )


def info(db):
    return dict(((r, k), v) for r, k, v in sqlite3.connect(db).execute("SELECT * FROM release_info"))


def main():
    schema = (ROOT / "assets/sql/db_schema.sqlite").read_text()
    releases = [(json.load(open(m))["release"], json.load(open(m))["date"]) for m in MANIFESTS]
    with tempfile.TemporaryDirectory() as tmp:
        db = f"{tmp}/a.sqlite3"
        make_db(db, schema, releases)
        out = run(db, *MANIFESTS)
        assert out.returncode == 0, out.stderr
        got = info(db)
        assert got[("1.0", "databases.pfam")] == "Pfam-A 38.0"
        assert got[("1.1", "databases.pfam")] == "Pfam-A 38.2"
        assert got[("1.1", "databases.pfam_profile_profile")] == "Pfam 35.0"
        assert got[("1.1", "databases.foldseek.pdb.pdb_date")] == "2025-01-01"
        assert got[("1.1", "n_families")] == "35459"
        assert not any(k.startswith("files") for _, k in got), "files list loaded"
        # Reload replaces, not duplicates
        assert run(db, MANIFESTS[1]).returncode == 0 and info(db) == got

        # A release missing from release_history, or a different date: refused, nothing changes
        db2 = f"{tmp}/b.sqlite3"
        make_db(db2, schema, releases[:1])
        out = run(db2, *MANIFESTS)
        assert out.returncode != 0 and "not in release_history" in out.stderr and info(db2) == {}, out.stderr
        db3 = f"{tmp}/c.sqlite3"
        make_db(db3, schema, [releases[0], (releases[1][0], "2026-11-17")])
        out = run(db3, *MANIFESTS)
        assert out.returncode != 0 and "date" in out.stderr and info(db3) == {}, out.stderr

        # A DB without the table (built before 3.1.0) gets it
        db4 = f"{tmp}/d.sqlite3"
        make_db(db4, schema.split("-- Metadata of each release")[0], releases)
        assert run(db4, *MANIFESTS).returncode == 0 and info(db4) == got
    print("test_load_release_manifest: OK")


if __name__ == "__main__":
    main()
