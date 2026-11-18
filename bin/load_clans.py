#!/usr/bin/env python3
"""Load a release's clans (FTP clans.tsv.gz, clan_edges.tsv.gz) into the mgnifam_clans and mgnifam_clan_edges tables
of an MGnifams sqlite DB, replacing the clans it held. One transaction; nothing changes unless every family is in the
DB, each clan has one representative, and the clans are exactly the connected components (two or more families) of
the edges with residue Jaccard >= --clan-jaccard. The tables and their indexes are created if missing.

    python3 bin/load_clans.py --db mgnifams.sqlite3 --clans clans.tsv.gz --edges clan_edges.tsv.gz
"""

import argparse
import csv
import gzip
import sqlite3
import sys
from collections import Counter

DDL = """
CREATE TABLE IF NOT EXISTS mgnifam_clans (
    mgnifam_id INTEGER PRIMARY KEY REFERENCES mgnifam(id),
    clan_id TEXT NOT NULL,
    is_representative INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS mgnifam_clan_edges (
    family_1 INTEGER NOT NULL REFERENCES mgnifam(id),
    family_2 INTEGER NOT NULL REFERENCES mgnifam(id),
    shared_proteins INTEGER,
    protein_jaccard REAL,
    residue_jaccard REAL,
    residue_containment REAL,
    PRIMARY KEY (family_1, family_2)
);
CREATE INDEX IF NOT EXISTS idx_clans_clan ON mgnifam_clans (clan_id);
CREATE INDEX IF NOT EXISTS idx_clan_edges_family_2 ON mgnifam_clan_edges (family_2);
"""
REPRESENTATIVE = {"True": 1, "False": 0}


def read_tsv(path, columns):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(columns) - set(reader.fieldnames or ())
        if missing:
            sys.exit(f"{path} lacks columns {sorted(missing)}")
        return [[row[c] for c in columns] for row in reader]


def components(edges, threshold):
    """Families grouped by the connected components of the edges with residue_jaccard >= threshold."""
    parent = {}

    def root(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for f1, f2, _, _, jaccard, _ in edges:
        if jaccard >= threshold:
            parent[root(f1)] = root(f2)
    groups = {}
    for family in parent:
        groups.setdefault(root(family), set()).add(family)
    return {frozenset(g) for g in groups.values() if len(g) > 1}


def check(clans, edges, threshold):
    reps = Counter(clan for clan, _, rep in clans if rep)
    members = {}
    for clan, family, _ in clans:
        members.setdefault(clan, set()).add(family)
    bad = [clan for clan in members if reps[clan] != 1]
    if bad:
        sys.exit(f"clans without exactly one representative: {bad[:20]}")
    if {frozenset(m) for m in members.values()} != components(edges, threshold):
        sys.exit(f"the clans are not the components of the edges with residue Jaccard >= {threshold}")


def load(con, clans, edges):
    con.execute("BEGIN IMMEDIATE")
    try:
        for statement in DDL.split(";"):
            if statement.strip():
                con.execute(statement)
        con.execute("DELETE FROM mgnifam_clans")
        con.execute("DELETE FROM mgnifam_clan_edges")
        con.executemany("INSERT INTO mgnifam_clans (clan_id, mgnifam_id, is_representative) VALUES (?, ?, ?)", clans)
        con.executemany("INSERT INTO mgnifam_clan_edges VALUES (?, ?, ?, ?, ?, ?)", edges)
        unknown = con.execute(
            "SELECT mgnifam_id FROM mgnifam_clans UNION SELECT family_1 FROM mgnifam_clan_edges "
            "UNION SELECT family_2 FROM mgnifam_clan_edges EXCEPT SELECT id FROM mgnifam LIMIT 20"
        ).fetchall()
        if unknown:
            raise ValueError(f"families missing from mgnifam: {[r[0] for r in unknown]}")
        con.execute("COMMIT")
    except (sqlite3.Error, ValueError) as e:
        con.execute("ROLLBACK")
        sys.exit(f"nothing loaded: {e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--clans", required=True, help="the release's clans.tsv.gz")
    ap.add_argument("--edges", required=True, help="the release's clan_edges.tsv.gz")
    ap.add_argument(
        "--clan-jaccard", type=float, default=0.5, help="the release's clan threshold (manifest clans.clan_jaccard)"
    )
    a = ap.parse_args()
    try:
        clans = [
            (clan, int(family), REPRESENTATIVE[rep])
            for clan, family, rep in read_tsv(a.clans, ("clan_id", "family_id", "is_representative"))
        ]
        edges = [
            (int(f1), int(f2), int(shared), float(pj), float(rj), float(rc))
            for f1, f2, shared, pj, rj, rc in read_tsv(
                a.edges,
                (
                    "family_1",
                    "family_2",
                    "shared_proteins",
                    "protein_jaccard",
                    "residue_jaccard",
                    "residue_containment",
                ),
            )
        ]
    except (KeyError, ValueError) as e:
        sys.exit(f"unreadable value: {e}")
    check(clans, edges, a.clan_jaccard)
    con = sqlite3.connect(a.db, isolation_level=None)
    load(con, clans, edges)
    n_clans = con.execute("SELECT count(DISTINCT clan_id) FROM mgnifam_clans").fetchone()[0]
    print(f"loaded {n_clans} clans of {len(clans)} families, {len(edges)} edges")


if __name__ == "__main__":
    main()
