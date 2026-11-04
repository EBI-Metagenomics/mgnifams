#!/usr/bin/env python3
"""Self-check for assets/merge_update_delta.sql (needs the sqlite3 CLI). Run: python3 bin/test_merge_update_delta.py"""

import hashlib
import sqlite3
import subprocess
import tempfile
from pathlib import Path

ASSETS = Path(__file__).resolve().parent.parent / "assets"
SCHEMA = (ASSETS / "data/db_schema.sqlite").read_text()


def make_db(path, rows, info=None):
    con = sqlite3.connect(path)
    con.executescript(SCHEMA)
    for i, full_size, blob in rows:
        con.execute(
            "INSERT INTO mgnifam (id, full_size, cif_blob, tm_blob, biome_blob, seed_msa_blob, seed_size, s4pred_blob,"
            " coil_percent) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (i, full_size, blob, blob, blob, blob if info is None else None, 5 if info is None else None, blob, full_size),
        )
        con.execute(
            "INSERT INTO mgnifam_pfams (mgnifam_id, pfam, e_value) VALUES (?, ?, 1e-5)", (i, f"PF{i}_{blob.decode()}")
        )
        con.execute(
            "INSERT INTO mgnifam_folds (mgnifam_id, fold, q_tmscore) VALUES (?, ?, 0.5)", (i, f"f_{blob.decode()}")
        )
    if info is not None:
        con.execute("CREATE TABLE update_info (key TEXT PRIMARY KEY, value TEXT)")
        con.executemany("INSERT INTO update_info VALUES (?, ?)", info.items())
    con.commit()
    con.close()


def merge(prod, delta):
    sql = (ASSETS / "merge_update_delta.sql").read_text()
    cmd = ["sqlite3", "-bail", prod, "-cmd", f"ATTACH '{delta}' AS delta"]
    return subprocess.run(cmd, input=sql, capture_output=True, text=True)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


INFO = {
    "tm_computed": "false",
    "biome_computed": "true",
    "child_tables": "pfams,funfams,folds",
    "mgnifams_release": "1.1",
    "mgnify_proteins_release": "2026_07",
    "pipeline_version": "3.0.0",
}

with tempfile.TemporaryDirectory() as d:
    prod, delta = f"{d}/prod.sqlite3", f"{d}/delta.sqlite3"
    make_db(prod, [(1, 10, b"old"), (2, 20, b"old"), (3, 30, b"old")])

    make_db(delta, [(1, 11, b"new"), (4, 44, b"new")], INFO)  # id 4 is not in prod
    before = sha(prod)
    r = merge(prod, delta)
    assert r.returncode != 0 and "delta_has_ids_missing_from_prod" in r.stderr, r.stderr
    assert sha(prod) == before, "failed merge must leave prod unchanged"

    Path(delta).unlink()
    make_db(delta, [(1, 11, b"new"), (2, 22, b"new")], {k: v for k, v in INFO.items() if k != "child_tables"})
    r = merge(prod, delta)
    assert r.returncode != 0 and "delta_lacks_update_info_child_tables" in r.stderr, r.stderr
    assert sha(prod) == before

    Path(delta).unlink()
    make_db(delta, [(1, 11, b"new"), (2, 22, b"new")], INFO)
    sqlite3.connect(delta).execute("DROP TABLE mgnifam_folds")
    r = merge(prod, delta)
    assert r.returncode != 0 and "no such table: delta.mgnifam_folds" in r.stderr, r.stderr
    assert sha(prod) == before, "a delta missing a child table must leave prod unchanged"

    Path(delta).unlink()
    make_db(delta, [(1, 11, b"new"), (2, 22, b"new")], INFO)
    with sqlite3.connect(delta) as con:
        con.execute("DELETE FROM mgnifam_folds WHERE mgnifam_id = 2")
    with sqlite3.connect(prod) as con:
        con.execute("UPDATE mgnifam SET has_pfam = 1, has_funfam = 1, has_model_pfam = 1, has_structure = 1")
    r = merge(prod, delta)
    assert r.returncode == 0, r.stderr
    con = sqlite3.connect(prod)
    rows = con.execute(
        "SELECT id, full_size, cif_blob, tm_blob, biome_blob, seed_msa_blob, seed_size FROM mgnifam ORDER BY id"
    ).fetchall()
    assert rows == [
        (1, 11, b"new", b"old", b"new", b"old", 5),  # tm_computed=false keeps tm_blob; seed MSA/size never change
        (2, 22, b"new", b"old", b"new", b"old", 5),
        (3, 30, b"old", b"old", b"old", b"old", 5),
    ], rows
    pfams = con.execute("SELECT mgnifam_id, pfam FROM mgnifam_pfams ORDER BY mgnifam_id").fetchall()
    assert pfams == [(1, "PF1_new"), (2, "PF2_new"), (3, "PF3_old")], pfams
    folds = con.execute("SELECT mgnifam_id, fold, q_tmscore FROM mgnifam_folds ORDER BY mgnifam_id").fetchall()
    assert folds == [(1, "f_new", 0.5), (3, "f_old", 0.5)], folds
    flags = con.execute(
        "SELECT id, has_pfam, has_funfam, has_model_pfam, has_structure FROM mgnifam ORDER BY id"
    ).fetchall()
    # Recomputed for the delta families only (no funfams, family 2 lost its fold); has_model_pfam is kept
    assert flags == [(1, 1, 0, 1, 1), (2, 1, 0, 1, 0), (3, 1, 1, 1, 1)], flags
    history = con.execute(
        "SELECT release, type, mgnify_proteins_release, pipeline_version, n_families, n_updated, n_not_updated FROM release_history"
    ).fetchall()
    assert history == [("1.1", "update", "2026_07", "3.0.0", 3, 2, 1)], history
    con.close()

    # The same release cannot be merged twice
    merged = sha(prod)
    r = merge(prod, delta)
    assert r.returncode != 0 and "UNIQUE constraint failed: release_history.release" in r.stderr, r.stderr
    assert sha(prod) == merged, "a repeated release must leave prod unchanged"

    # A delta without the release keys is rejected; a prod without release_history gets it created
    for key in ("mgnifams_release", "mgnify_proteins_release"):
        Path(delta).unlink()
        make_db(delta, [(1, 11, b"new")], {k: v for k, v in INFO.items() if k != key})
        r = merge(prod, delta)
        assert r.returncode != 0 and "NOT NULL constraint failed: release_history" in r.stderr, r.stderr
        assert sha(prod) == merged

    Path(prod).unlink()
    make_db(prod, [(1, 10, b"old")])
    with sqlite3.connect(prod) as con:
        con.execute("DROP TABLE release_history")
    Path(delta).unlink()
    make_db(delta, [(1, 11, b"new")], {**INFO, "mgnifams_release": "1.2"})
    r = merge(prod, delta)
    assert r.returncode == 0, r.stderr
    con = sqlite3.connect(prod)
    assert con.execute("SELECT release, n_not_updated FROM release_history").fetchall() == [("1.2", 0)]
    con.close()

# Annotation-only deltas (annotate_not_updated): only s4pred, TM (when computed) and the child tables change
INFO_ANN = {**INFO, "annotation_only": "true", "tm_computed": "false", "afdb_computed": "true"}
KEPT = "SELECT id, full_size, cif_blob, tm_blob, biome_blob, seed_msa_blob, seed_size FROM mgnifam ORDER BY id"


def content(path):
    """Rows that matter, without the child tables' AUTOINCREMENT ids (re-inserted rows get new ones)."""
    con = sqlite3.connect(path)
    out = [
        con.execute("SELECT * FROM mgnifam ORDER BY id").fetchall(),
        con.execute("SELECT mgnifam_id, pfam, e_value FROM mgnifam_pfams ORDER BY mgnifam_id, pfam").fetchall(),
        con.execute("SELECT mgnifam_id, fold, q_tmscore FROM mgnifam_folds ORDER BY mgnifam_id, fold").fetchall(),
        con.execute("SELECT * FROM release_history ORDER BY release").fetchall(),
    ]
    con.close()
    return out


with tempfile.TemporaryDirectory() as d:
    prod, delta, ann = f"{d}/prod.sqlite3", f"{d}/delta.sqlite3", f"{d}/ann.sqlite3"
    make_db(prod, [(1, 10, b"old"), (2, 20, b"old"), (3, 30, b"old")])
    make_db(delta, [(1, 11, b"upd")], INFO)  # the release's update delta goes first
    assert merge(prod, delta).returncode == 0
    con = sqlite3.connect(prod)
    kept_before = con.execute(KEPT).fetchall()
    history_before = con.execute("SELECT * FROM release_history").fetchall()
    con.close()

    make_db(ann, [(2, 22, b"new"), (3, 33, b"new")], INFO_ANN)
    r = merge(prod, ann)
    assert r.returncode == 0, r.stderr
    con = sqlite3.connect(prod)
    # ann_keeps_columns: rep / size / structure / TM (not computed) / biome / seed columns untouched
    assert con.execute(KEPT).fetchall() == kept_before, con.execute(KEPT).fetchall()
    # ann_updates_s4pred
    s4 = con.execute("SELECT id, s4pred_blob, coil_percent FROM mgnifam ORDER BY id").fetchall()
    assert s4 == [(1, b"upd", 11.0), (2, b"new", 22.0), (3, b"new", 33.0)], s4
    # ann_child_tables: replaced for 2 and 3, family 1 keeps the update's rows; flags recomputed
    pfams = con.execute("SELECT mgnifam_id, pfam FROM mgnifam_pfams ORDER BY mgnifam_id").fetchall()
    assert pfams == [(1, "PF1_upd"), (2, "PF2_new"), (3, "PF3_new")], pfams
    flags = con.execute("SELECT id, has_pfam, has_funfam, has_structure FROM mgnifam ORDER BY id").fetchall()
    assert flags == [(1, 1, 0, 1), (2, 1, 0, 1), (3, 1, 0, 1)], flags
    # ann_no_history
    assert con.execute("SELECT * FROM release_history").fetchall() == history_before
    con.close()

    # ann_twice: a repeated annotation-only merge succeeds and changes no content
    once = content(prod)
    r = merge(prod, ann)
    assert r.returncode == 0, r.stderr
    assert content(prod) == once

    # An annotation-only delta on a prod without any update merged: no release_history row
    Path(prod).unlink()
    make_db(prod, [(2, 20, b"old"), (3, 30, b"old")])
    r = merge(prod, ann)
    assert r.returncode == 0, r.stderr
    assert sqlite3.connect(prod).execute("SELECT count(*) FROM release_history").fetchone()[0] == 0

# ann_afdb_required: an annotation-only delta without the AlphaFoldDB search would drop prod's AFDB folds
with tempfile.TemporaryDirectory() as d:
    prod, ann = f"{d}/prod.sqlite3", f"{d}/ann.sqlite3"
    make_db(prod, [(2, 20, b"old"), (3, 30, b"old")])
    before = sha(prod)
    for afdb in ("false", None):
        Path(ann).unlink(missing_ok=True)
        info = {k: v for k, v in INFO_ANN.items() if k != "afdb_computed"} | ({"afdb_computed": afdb} if afdb else {})
        make_db(ann, [(2, 22, b"new")], info)
        r = merge(prod, ann)
        assert r.returncode != 0 and "annotation_delta_lacks_afdb" in r.stderr, (afdb, r.stderr)
        assert sha(prod) == before, "a rejected annotation-only delta must leave prod unchanged"
print("ok")
