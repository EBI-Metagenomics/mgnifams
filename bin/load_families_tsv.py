#!/usr/bin/env python3
"""Load a release's families.tsv.gz provenance columns into the mgnifam table of an MGnifams sqlite DB.

Run after the release's merges. One transaction: every family of the DB gets the values of its row, and the DB and
the file must hold the same family ids, or nothing changes. An empty not_updated_reason is stored as NULL.

    python3 bin/load_families_tsv.py --db mgnifams.sqlite3 --families families.tsv.gz
"""

import argparse
import csv
import gzip
import sqlite3
import sys

COLUMNS = (
    "status",
    "not_updated_reason",
    "first_release",
    "model_release",
    "model_proteins_release",
    "members_release",
    "members_proteins_release",
)


def read_rows(path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"family_id", *COLUMNS} - set(reader.fieldnames or ())
        if missing:
            sys.exit(f"{path} lacks columns {sorted(missing)}")
        rows = [(int(row["family_id"]), *(row[c] or None for c in COLUMNS)) for row in reader]
    empty = [r[0] for r in rows if any(v is None for c, v in zip(COLUMNS, r[1:]) if c != "not_updated_reason")]
    if empty:
        sys.exit(f"{path}: families with empty provenance values: {empty[:20]}")
    return rows


def load(con, rows):
    have = {r[1] for r in con.execute("SELECT * FROM pragma_table_info('mgnifam')")}
    if set(COLUMNS) - have:
        sys.exit("mgnifam lacks the provenance columns: run assets/sql/migrate_schema_family_status.sql first")
    con.execute("BEGIN IMMEDIATE")
    try:
        con.execute(f"CREATE TEMP TABLE families (id INTEGER PRIMARY KEY, {', '.join(COLUMNS)})")
        con.executemany(f"INSERT INTO temp.families VALUES ({', '.join('?' * (len(COLUMNS) + 1))})", rows)
        only_db = con.execute("SELECT id FROM main.mgnifam EXCEPT SELECT id FROM temp.families LIMIT 20").fetchall()
        only_tsv = con.execute("SELECT id FROM temp.families EXCEPT SELECT id FROM main.mgnifam LIMIT 20").fetchall()
        if only_db or only_tsv:
            raise ValueError(f"family ids differ: only in the DB {only_db}, only in the file {only_tsv}")
        con.execute(
            f"UPDATE main.mgnifam SET {', '.join(f'{c} = f.{c}' for c in COLUMNS)} "
            "FROM temp.families AS f WHERE main.mgnifam.id = f.id"
        )
        con.execute("DROP TABLE temp.families")
        con.execute("COMMIT")
    except (sqlite3.Error, ValueError) as e:
        con.execute("ROLLBACK")
        sys.exit(f"nothing loaded: {e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--families", required=True, help="the release's families.tsv.gz")
    a = ap.parse_args()
    rows = read_rows(a.families)
    con = sqlite3.connect(a.db, isolation_level=None)
    load(con, rows)
    counts = con.execute("SELECT status, count(*) FROM mgnifam GROUP BY status ORDER BY status").fetchall()
    print(f"loaded {len(rows)} families: {', '.join(f'{s} {n}' for s, n in counts)}")


if __name__ == "__main__":
    main()
