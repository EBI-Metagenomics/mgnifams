#!/usr/bin/env python3
"""Load per-family blob files into an update (delta) MGnifams sqlite DB, in one transaction, and refuse to leave
a partial DB behind: exits non-zero (and rolls back) on a duplicate or empty file, a missing required blob, an empty
always-updated scalar, or a family id set that differs from the expected one.

Files are found recursively under each directory as <family_id>.<ext>.
"""

import argparse
import sqlite3
import sys
from pathlib import Path

BLOBS = {  # name: (column, extension)
    "cif": ("cif_blob", ".cif"),
    "s4pred": ("s4pred_blob", ".json"),
    "tm": ("tm_blob", ".json"),
    "biome": ("biome_blob", ".csv"),
    "domain": ("domain_blob", ".json"),
}
SCALARS = "full_size protein_rep rep_region rep_length rep_sequence plddt ptm helix_percent strand_percent coil_percent"
TM_SCALARS = (
    "inside_percent membrane_alpha_percent outside_percent signal_percent membrane_beta_percent periplasm_percent"
)


def fail(msg):
    raise SystemExit(f"update_sqlite_blobs_staged: {msg}")


def blob_files(directory, ext):
    files = {}
    for f in sorted(Path(directory).rglob(f"*{ext}")):
        family = f.name[: -len(ext)]
        if family in files:
            fail(f"duplicate {f.name}: {files[family]} and {f}")
        files[family] = f
    return files


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True)
    for name in BLOBS:
        parser.add_argument(f"--{name}_dir", help=f"Directory with <id>{BLOBS[name][1]} files")
    parser.add_argument("--required", default="", help="Comma-separated blobs every family must have, e.g. cif,s4pred")
    parser.add_argument("--ids", required=True, help="File with the expected family ids, one per line")
    parser.add_argument("--info", nargs="*", default=[], help="key=value rows for the update_info table")
    a = parser.parse_args(args)

    required = [r for r in a.required.split(",") if r]
    expected = {line.strip() for line in open(a.ids) if line.strip()}
    con = sqlite3.connect(a.db)
    try:
        with con:  # one transaction; any exception rolls it back
            ids = {str(r[0]) for r in con.execute("SELECT id FROM mgnifam")}
            if ids != expected:
                fail(
                    f"mgnifam ids differ from the expected ids: missing {sorted(expected - ids)[:20]}, unexpected {sorted(ids - expected)[:20]}"
                )
            for name, (column, ext) in BLOBS.items():
                directory = getattr(a, f"{name}_dir")
                if not directory:
                    continue
                for family, f in blob_files(directory, ext).items():
                    if family not in ids:
                        continue
                    data = f.read_bytes()
                    if not data:
                        fail(f"empty file {f}")
                    con.execute(f"UPDATE mgnifam SET {column} = ? WHERE id = ?", (data, family))

            checks = [f"{BLOBS[r][0]} IS NULL OR length({BLOBS[r][0]}) = 0" for r in required]
            scalars = SCALARS.split() + (TM_SCALARS.split() if "tm" in required else [])
            checks += [f"{c} IS NULL OR {c} = ''" for c in scalars]
            bad = [str(r[0]) for r in con.execute(f"SELECT id FROM mgnifam WHERE {' OR '.join(checks)} ORDER BY id")]
            if bad:
                fail(
                    f"{len(bad)} families lack a required blob or value ({', '.join(required + ['scalars'])}): {bad[:20]}"
                )

            con.execute("CREATE TABLE IF NOT EXISTS update_info (key TEXT PRIMARY KEY, value TEXT)")
            con.executemany("INSERT OR REPLACE INTO update_info VALUES (?, ?)", [kv.split("=", 1) for kv in a.info])
    except SystemExit as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    finally:
        con.close()


if __name__ == "__main__":
    main()
