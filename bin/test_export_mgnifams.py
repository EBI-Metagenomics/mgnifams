#!/usr/bin/env python3
"""Self-check: mgnifam.csv from export_mgnifams.py lands in the right columns when IMPORT_QUERIES imports it
by position into assets/data/db_schema.sqlite. Needs the sqlite3 CLI. Run: python3 bin/test_export_mgnifams.py"""

import sqlite3
import subprocess
import tempfile
from pathlib import Path

from export_mgnifams import write_mgnifam_csv

SCHEMA = (Path(__file__).resolve().parent.parent / "assets/data/db_schema.sqlite").read_text()

with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    (d / "meta.csv").write_text(
        "family_id,full_msa_size,protein,region,length,sequence,consensus,converged\n"
        '1,21,"33",357-470,114,QLDP,lldp,True\n'
    )
    (d / "scores.csv").write_text("id,plddt,ptm\n1,80.5,0.7\n")
    (d / "comp.csv").write_text("id,helix_percent,strand_percent,coil_percent\n1,10.0,20.0,70.0\n")
    (d / "seeds.csv").write_text("id,seed_size\n1,14\n")
    write_mgnifam_csv(d / "meta.csv", d / "scores.csv", d / "comp.csv", "", d / "mgnifam.csv", d / "seeds.csv")

    db = d / "db.sqlite3"
    sqlite3.connect(db).executescript(SCHEMA)
    subprocess.run(["sqlite3", db, ".import --csv --skip 1 mgnifam.csv mgnifam"], cwd=d, check=True)
    row = (
        sqlite3.connect(db)
        .execute(
            "SELECT full_size, seed_size, converged, protein_rep, rep_region, rep_length, rep_sequence, consensus, hmm_length,"
            " plddt, ptm, coil_percent FROM mgnifam WHERE id = 1"
        )
        .fetchone()
    )
    assert row == (21, 14, "True", 33, "357-470", 114, "QLDP", "lldp", 4, 80.5, 0.7, 70.0), row

print("ok")
