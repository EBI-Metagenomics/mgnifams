#!/usr/bin/env python3
"""Self-check for bin/load_clans.py. Run: python3 bin/test_load_clans.py"""

import gzip
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = Path(__file__).resolve().parent
SCHEMA = (BIN.parent / "assets/sql/db_schema.sqlite").read_text()
CLANS = "clan_id\tfamily_id\tis_representative\n"
EDGES = "family_1\tfamily_2\tshared_proteins\tprotein_jaccard\tresidue_jaccard\tresidue_containment\n"
# Families 1-2-3 joined at Jaccard >= 0.5 (one clan), 3-4 only a weak edge, 5 alone
GOOD_CLANS = ["1.1_clan_1\t1\tTrue\n", "1.1_clan_1\t2\tFalse\n", "1.1_clan_1\t3\tFalse\n"]
GOOD_EDGES = ["1\t2\t10\t0.9\t0.8\t0.9\n", "2\t3\t5\t0.6\t0.5\t0.7\n", "3\t4\t2\t0.1\t0.2\t0.6\n"]


def write(path, header, rows):
    with gzip.open(path, "wt") as handle:
        handle.write(header + "".join(rows))


def run(tmp, db, clans, edges):
    write(f"{tmp}/c.tsv.gz", CLANS, clans)
    write(f"{tmp}/e.tsv.gz", EDGES, edges)
    return subprocess.run(
        [
            sys.executable,
            str(BIN / "load_clans.py"),
            "--db",
            db,
            "--clans",
            f"{tmp}/c.tsv.gz",
            "--edges",
            f"{tmp}/e.tsv.gz",
        ],
        capture_output=True,
        text=True,
    )


def state(db):
    con = sqlite3.connect(db)
    return (
        con.execute("SELECT * FROM mgnifam_clans ORDER BY mgnifam_id").fetchall(),
        con.execute("SELECT count(*) FROM mgnifam_clan_edges").fetchone()[0],
    )


def main():
    with tempfile.TemporaryDirectory() as tmp:
        db = f"{tmp}/a.sqlite3"
        con = sqlite3.connect(db)
        con.executescript(SCHEMA)
        con.executemany("INSERT INTO mgnifam (id) VALUES (?)", [(i,) for i in range(1, 6)])
        con.commit()
        con.close()
        out = run(tmp, db, GOOD_CLANS, GOOD_EDGES)
        assert out.returncode == 0 and "loaded 1 clans of 3 families, 3 edges" in out.stdout, out.stderr + out.stdout
        want = ([(1, "1.1_clan_1", 1), (2, "1.1_clan_1", 0), (3, "1.1_clan_1", 0)], 3)
        assert state(db) == want, state(db)
        # Reload replaces
        assert run(tmp, db, GOOD_CLANS, GOOD_EDGES).returncode == 0 and state(db) == want

        failing = {
            "two representatives": (GOOD_CLANS[:1] + ["1.1_clan_1\t2\tTrue\n"] + GOOD_CLANS[2:], GOOD_EDGES),
            "clan is not a component": (GOOD_CLANS + ["1.1_clan_1\t4\tFalse\n"], GOOD_EDGES),
            "component missing": (GOOD_CLANS[:2], GOOD_EDGES),
            "unknown family": (GOOD_CLANS, GOOD_EDGES + ["3\t9\t1\t0.1\t0.1\t0.5\n"]),
            "duplicate edge": (GOOD_CLANS, GOOD_EDGES + GOOD_EDGES[2:]),
        }
        for name, (clans, edges) in failing.items():
            out = run(tmp, db, clans, edges)
            assert out.returncode != 0, name
            assert state(db) == want, name

        # A DB without the tables gets them
        db2 = f"{tmp}/b.sqlite3"
        con = sqlite3.connect(db2)
        con.executescript(SCHEMA.split("-- Clans of the release")[0])
        con.executemany("INSERT INTO mgnifam (id) VALUES (?)", [(i,) for i in range(1, 6)])
        con.commit()
        con.close()
        assert run(tmp, db2, GOOD_CLANS, GOOD_EDGES).returncode == 0 and state(db2) == want
    print("test_load_clans: OK")


if __name__ == "__main__":
    main()
