#!/usr/bin/env python3
"""Self-check: mgnifam.csv from export_mgnifams.py imports into assets/data/db_schema.sqlite the way IMPORT_QUERIES
does it (by header name), so every header column exists in the schema. Needs the sqlite3 CLI.
Run: python3 bin/test_export_mgnifams.py"""

import sqlite3
import subprocess
import tempfile
from pathlib import Path

from export_mgnifams import write_mgnifam_csv

SCHEMA = (Path(__file__).resolve().parent.parent / "assets/data/db_schema.sqlite").read_text()

QUERY = (
    "SELECT full_size, seed_size, converged, protein_rep, rep_region, rep_length, rep_sequence, consensus, hmm_length,"
    " plddt, ptm, coil_percent, has_pfam, cif_blob FROM mgnifam WHERE id = 1"
)


def imported_row(metadata, seed_sizes):
    """Export one family and import it; `seed_sizes` is the id,seed_size CSV text, or None."""
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "meta.csv").write_text(metadata)
        (d / "scores.csv").write_text("id,plddt,ptm\n1,80.5,0.7\n")
        (d / "comp.csv").write_text("id,helix_percent,strand_percent,coil_percent\n1,10.0,20.0,70.0\n")
        seeds = ""
        if seed_sizes is not None:
            seeds = d / "seeds.csv"
            seeds.write_text(seed_sizes)
        write_mgnifam_csv(d / "meta.csv", d / "scores.csv", d / "comp.csv", "", d / "mgnifam.csv", seeds)

        db = d / "db.sqlite3"
        sqlite3.connect(db).executescript(SCHEMA)
        cols = (d / "mgnifam.csv").read_text().splitlines()[0]
        nullif = ", ".join(f"NULLIF({c}, '')" for c in cols.split(","))
        sql = f".bail on\n.import --csv mgnifam.csv temp_mgnifam\nINSERT INTO mgnifam ({cols}) SELECT {nullif} FROM temp_mgnifam;\n"
        subprocess.run(["sqlite3", db], input=sql, text=True, cwd=d, check=True)
        return sqlite3.connect(db).execute(QUERY).fetchone()


# Default mode: bin/generate_families.py metadata, seed sizes from COUNT_SEED_MSA_SIZES, hmm_length from the consensus.
# has_* keep their DEFAULT (FINALIZE_SQLITE fills them); empty blob cells become NULL
row = imported_row(
    "family_id,full_msa_size,protein,region,length,sequence,consensus,converged\n"
    '1,21,"33",357-470,114,QLDP,lldp,True\n',
    "id,seed_size\n1,14\n",
)
assert row == (21, 14, "True", 33, "357-470", 114, "QLDP", "lldp", 4, 80.5, 0.7, 70.0, 0, None), row

# Update mode: mgnifam update_families --skip_refine metadata, read by name; seed size and converged are empty
row = imported_row(
    "family_id,converged,seed_msa_size,full_msa_size,rep_protein,rep_region,rep_length,consensus_length,rep_sequence,consensus_sequence\n"
    '1,,,21,"33",357-470,114,4,QLDP,lldp\n',
    None,
)
assert row == (21, None, None, 33, "357-470", 114, "QLDP", "lldp", 4, 80.5, 0.7, 70.0, 0, None), row

print("ok")
