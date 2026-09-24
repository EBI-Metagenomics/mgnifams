#!/usr/bin/env python3
"""Write an update release's families.tsv.gz (columns in ftp/README.md) from the previous release's file.

Every previous family is kept: `updated` rows take the new members from family_metadata.csv, `not_updated`
rows keep their previous values with the update's discard reason (updated_delta.csv `outcome`).
"""

import argparse
import csv
import gzip
import io
import sys


def read_delta(path):
    # mgnifam writes the delta unquoted, with `outcome` last: split on the first 8 commas only
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split(",")
        return {
            fields[0]: fields[-1]
            for fields in (line.rstrip("\n").split(",", len(header) - 1) for line in handle if line.strip())
        }


def update_rows(previous, outcomes, metadata, release, proteins_release):
    ids = [row["family_id"] for row in previous]
    for label, found in (("updated_delta.csv", outcomes.keys()), ("family_metadata.csv", metadata.keys())):
        unknown = set(found) - set(ids)
        if unknown:
            sys.exit(f"{label} has families missing from the previous families.tsv: {sorted(unknown)[:20]}")
    missing = set(ids) - set(outcomes)
    if missing:
        sys.exit(f"previous families missing from updated_delta.csv: {sorted(missing)[:20]}")
    successful = {family for family, outcome in outcomes.items() if outcome == "successful"}
    if successful != set(metadata):
        sys.exit(f"successful families without metadata: {sorted(successful ^ set(metadata))[:20]}")
    if any(row["members_release"] == release for row in previous):
        sys.exit(f"the previous families.tsv already has members from release {release}")

    rows = []
    for row in previous:
        family = row["family_id"]
        if family in successful:
            meta = metadata[family]
            row = row | {
                "status": "updated",
                "not_updated_reason": "",
                "members_release": release,
                "members_proteins_release": proteins_release,
                "full_size": meta["full_msa_size"],
                "rep_id": f"MGYP{int(meta['protein']):012d}",
                "rep_region": meta["region"],
            }
        else:
            row = row | {"status": "not_updated", "not_updated_reason": outcomes[family]}
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, help="previous release's families.tsv.gz")
    parser.add_argument("--delta", required=True, help="updated_delta.csv")
    parser.add_argument("--metadata", required=True, help="family_metadata.csv")
    parser.add_argument("--release", required=True, help="this MGnifams release, MAJOR.MINOR")
    parser.add_argument("--proteins_release", required=True, help="MGnify Proteins release of the members, YYYY_MM")
    parser.add_argument("--output", required=True, help="families.tsv.gz")
    args = parser.parse_args()

    with gzip.open(args.previous, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        previous = list(reader)
    with open(args.metadata, newline="") as handle:
        metadata = {row["family_id"]: row for row in csv.DictReader(handle)}

    rows = update_rows(previous, read_delta(args.delta), metadata, args.release, args.proteins_release)
    # mtime=0: the same inputs give a byte-identical file
    with gzip.GzipFile(args.output, "wb", mtime=0) as raw, io.TextIOWrapper(raw, newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
