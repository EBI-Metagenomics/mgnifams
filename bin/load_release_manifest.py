#!/usr/bin/env python3
"""Load FTP release manifests (ftp/releases/X.Y/release_manifest.json) into the release_info table of an MGnifams
sqlite DB, flattened to dotted keys (databases.pfam, databases.foldseek.pdb.pdb_date, ...); the per-file list is left
out. Each manifest replaces its release's rows. All manifests load in one transaction, and only for releases in
release_history whose date matches, or nothing changes. The table is created if missing.

    python3 bin/load_release_manifest.py --db mgnifams.sqlite3 --manifests releases/1.0/release_manifest.json releases/1.1/release_manifest.json
"""

import argparse
import json
import sqlite3
import sys

SKIP = {"files"}

DDL = """CREATE TABLE IF NOT EXISTS release_info (
    release TEXT NOT NULL REFERENCES release_history(release),
    key TEXT NOT NULL,
    value TEXT,
    PRIMARY KEY (release, key)
)"""


def flatten(node, prefix=""):
    for key, value in node.items():
        if not prefix and key in SKIP:
            continue
        if isinstance(value, dict):
            yield from flatten(value, f"{prefix}{key}.")
        else:
            yield f"{prefix}{key}", value if isinstance(value, str) else json.dumps(value)


def load(con, manifests):
    con.execute("BEGIN IMMEDIATE")
    try:
        con.execute(DDL)
        for path, manifest in manifests:
            release, date = manifest.get("release"), manifest.get("date")
            row = con.execute("SELECT date FROM release_history WHERE release = ?", (release,)).fetchone()
            if row is None:
                raise ValueError(f"{path}: release {release!r} is not in release_history")
            if row[0] != date:
                raise ValueError(f"{path}: date {date} != release_history date {row[0]} of release {release}")
            con.execute("DELETE FROM release_info WHERE release = ?", (release,))
            con.executemany("INSERT INTO release_info VALUES (?, ?, ?)", [(release, k, v) for k, v in flatten(manifest)])
        con.execute("COMMIT")
    except (sqlite3.Error, ValueError) as e:
        con.execute("ROLLBACK")
        sys.exit(f"nothing loaded: {e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--manifests", nargs="+", required=True)
    a = ap.parse_args()
    manifests = []
    for path in a.manifests:
        with open(path) as handle:
            manifests.append((path, json.load(handle)))
    con = sqlite3.connect(a.db, isolation_level=None)
    load(con, manifests)
    for release, n in con.execute("SELECT release, count(*) FROM release_info GROUP BY release ORDER BY release"):
        print(f"release_info {release}: {n} keys")


if __name__ == "__main__":
    main()
