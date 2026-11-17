#!/usr/bin/env python3
"""Self-check for bin/load_families_tsv.py and assets/sql/migrate_schema_family_status.sql.
Run: python3 bin/test_load_families_tsv.py"""

import gzip
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = Path(__file__).resolve().parent
SQL = BIN.parent / "assets/sql"
HEADER = "family_id\tstatus\tnot_updated_reason\tfirst_release\tmodel_release\tmodel_proteins_release\tmembers_release\tmembers_proteins_release\tfull_size\tseed_size\trep_id\trep_region\n"
ROWS = [
    "1\tupdated\t\t1.0\t1.0\t2024_04\t1.1\t2026_07\t10\t2\tMGYP000000000001\t1-50\n",
    "2\tnot_updated\tfamily representative length too small\t1.0\t1.0\t2024_04\t1.0\t2024_04\t5\t2\tMGYP000000000002\t-\n",
]


def write_tsv(path, rows):
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER + "".join(rows))


def make_db(path, schema):
    con = sqlite3.connect(path)
    con.executescript(schema)
    con.executemany("INSERT INTO mgnifam (id, full_size) VALUES (?, 1)", [(1,), (2,)])
    con.commit()
    con.close()


def run(db, tsv):
    return subprocess.run(
        [sys.executable, str(BIN / "load_families_tsv.py"), "--db", db, "--families", tsv],
        capture_output=True,
        text=True,
    )


def status(db):
    return sqlite3.connect(db).execute(
        "SELECT id, status, not_updated_reason, members_release FROM mgnifam ORDER BY id"
    ).fetchall()


def main():
    with tempfile.TemporaryDirectory() as tmp:
        db, tsv, bad = f"{tmp}/a.sqlite3", f"{tmp}/f.tsv.gz", f"{tmp}/bad.tsv.gz"
        make_db(db, (SQL / "db_schema.sqlite").read_text())
        write_tsv(tsv, ROWS)
        out = run(db, tsv)
        assert out.returncode == 0, out.stderr
        want = [(1, "updated", None, "1.1"), (2, "not_updated", "family representative length too small", "1.0")]
        assert status(db) == want, status(db)
        assert "not_updated 1, updated 1" in out.stdout, out.stdout

        # Ids differ (a family missing, an extra one): fails, the DB keeps the first load
        write_tsv(bad, [ROWS[0], ROWS[1].replace("2\t", "3\t", 1)])
        out = run(db, bad)
        assert out.returncode != 0 and "family ids differ" in out.stderr, out.stderr
        assert status(db) == want
        # Duplicate id: fails
        write_tsv(bad, ROWS + ROWS[:1])
        assert run(db, bad).returncode != 0 and status(db) == want

        # A DB without the columns: refused, then migrated (appended columns) and loaded
        old = (SQL / "db_schema.sqlite").read_text()
        old = "\n".join(line for line in old.splitlines() if not any(
            line.strip().startswith(c + " ") for c in ("status", "not_updated_reason", "first_release", "model_release",
                                                      "model_proteins_release", "members_release",
                                                      "members_proteins_release")))
        db2 = f"{tmp}/b.sqlite3"
        make_db(db2, old)
        out = run(db2, tsv)
        assert out.returncode != 0 and "migrate_schema_family_status.sql" in out.stderr, out.stderr
        con = sqlite3.connect(db2)
        con.executescript((SQL / "migrate_schema_family_status.sql").read_text())
        con.close()
        assert run(db2, tsv).returncode == 0 and status(db2) == want
    print("test_load_families_tsv: OK")


if __name__ == "__main__":
    main()
