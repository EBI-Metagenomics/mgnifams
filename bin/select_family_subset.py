#!/usr/bin/env python3
"""Pick the families of one status (e.g. not_updated) of a release's families.tsv.gz out of a reps FASTA and a
folder of <family_id>.pdb structures. FASTA records are matched on the first header token and written with the bare
family id; PDBs are symlinked. Exits 1, naming the ids, when a selected family is missing from either input; an
empty selection gives empty outputs (exit 0).
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path


def read_fasta(path):
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                records[name] = []
            elif name is not None and line:
                records[name].append(line)
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", required=True, help="families.tsv.gz of the release")
    parser.add_argument("--reps", required=True, help="FASTA whose headers start with the family id")
    parser.add_argument("--structures", required=True, help="folder of <family_id>.pdb")
    parser.add_argument("--status", required=True, help="families.tsv status to select, e.g. not_updated")
    parser.add_argument("--out_fasta", required=True)
    parser.add_argument("--out_structures", required=True)
    a = parser.parse_args(argv)

    with gzip.open(a.families, "rt", newline="") as fh:
        ids = [row["family_id"] for row in csv.DictReader(fh, delimiter="\t") if row["status"] == a.status]

    reps = read_fasta(a.reps)
    pdbs = {p.name[: -len(".pdb")]: p for p in Path(a.structures).iterdir() if p.name.endswith(".pdb")}
    for what, found in (("reps", reps), ("structures", pdbs)):
        missing = [i for i in ids if i not in found]
        if missing:
            print(f"select_family_subset: {len(missing)} families missing from {what}: {missing[:20]}", file=sys.stderr)
            return 1

    out = Path(a.out_structures)
    out.mkdir(parents=True, exist_ok=True)
    with open(a.out_fasta, "w") as fh:
        for i in ids:
            fh.write(f">{i}\n" + "".join(f"{s}\n" for s in reps[i]))
            (out / f"{i}.pdb").symlink_to(pdbs[i].resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
