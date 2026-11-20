#!/usr/bin/env python3
"""Classify the proteins of MGnifams full MSAs by how many copies of the family domain they carry.

Each full-MSA row is one HMMER domain region `MGYP/start-end`, aligned by hmmalign to the family HMM.
Its model span is the first and last match column (`#=GC RF` not ".") holding a residue.
Within one protein, rows in protein order are chained into copies: a row continues the previous
copy when its model span starts later and overlaps the previous row by at most 50% of the shorter
row (one domain split by an insertion, as InterPro's discontinuous domains); otherwise it starts
a new copy. A copy is full length when it spans >= 90% of the model from its first to its last model
position, as a single row is judged (the generation length filter).

Protein classes: single (one row), split (one copy from >= 2 rows), tandem (>= 2 copies, all full),
mixed (full and partial copies), partial_repeat (>= 2 copies, none full).
Row ends follow Pfam's graphics: curved when the match reaches within 5% of the model of position 1 / the last
position (strict 1 / M left every MGnifam 25 end jagged, its rows span 3-275 of 278).

Usage: classify_multi_hits.py MSA_OR_DIR... [--rows rows.tsv.gz] > summary.tsv
Reads `<family>.sto.gz` single-block Stockholm (mgnifam writes format="pfam"). Read-only.
"""

import argparse
import gzip
from pathlib import Path

FULL = 0.9
OVERLAP = 0.5
END = 0.05  # per end: the two ends together allow the 10% of the full-length rule
CLASSES = ("split", "tandem", "mixed", "partial_repeat")
FIELDS = (
    "family",
    "model_length",
    "regions",
    "mgyps",
    "full_length_regions",
    "single_full",
    "single_partial",
    "multi_hit_mgyps",
    "split",
    "split_full",
    "tandem",
    "mixed",
    "partial_repeat",
    "copies",
    "full_copies",
)
ROW_FIELDS = (
    "family",
    "mgyp",
    "seq_from",
    "seq_to",
    "hmm_from",
    "hmm_to",
    "n_end",
    "c_end",
    "copy",
    "copy_full",
    "class",
)


def model_span(aligned, rf):
    """First and last model position with a residue in this row, or None."""
    position, first, last = 0, None, None
    for residue, ref in zip(aligned, rf):
        if ref == ".":
            continue
        position += 1
        if residue.isalpha():
            first = first or position
            last = position
    return None if first is None else (first, last)


def region(name, aligned):
    """(mgyp, start, end) from `mgyp/start-end`; a name without "/" is the whole protein."""
    if "/" not in name:
        return name, 1, sum(c.isalpha() for c in aligned)
    mgyp, coords = name.rsplit("/", 1)
    start, end = coords.split("-")
    return mgyp, int(start), int(end)


def ends(hmm_from, hmm_to, model_length):
    """Pfam's graphics convention, with up to 5% of the model missing at an end still drawn curved."""
    tolerance = END * model_length
    return (
        "curved" if hmm_from - 1 <= tolerance else "jagged",
        "curved" if model_length - hmm_to <= tolerance else "jagged",
    )


def chain(rows, model_length, overlap=OVERLAP):
    """Group (seq_from, seq_to, hmm_from, hmm_to) rows of one protein into copies, in protein order."""
    copies = []
    for row in sorted(rows):
        if copies:
            last = copies[-1][-1]
            shared = last[3] - row[2] + 1
            if row[2] > last[2] and shared <= overlap * min(last[3] - last[2] + 1, row[3] - row[2] + 1):
                copies[-1].append(row)
                continue
        copies.append([row])
    return copies


def is_full(copy, model_length):
    """Judged by its ends, as a single row: positions skipped inside a split copy still count."""
    return max(row[3] for row in copy) - min(row[2] for row in copy) + 1 >= FULL * model_length


def classify(copies, model_length):
    if len(copies) == 1:
        return "single" if len(copies[0]) == 1 else "split"
    full = sum(is_full(copy, model_length) for copy in copies)
    return "tandem" if full == len(copies) else "mixed" if full else "partial_repeat"


def read_msa(path):
    """(family, rf, [(name, aligned)]) from a single-block Stockholm file."""
    rf, rows, seen = None, [], set()
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        for line in handle:
            if line.startswith("#=GC RF"):
                if rf is not None:
                    raise ValueError(f"{path}: more than one RF line (multi-block Stockholm)")
                rf = line.split()[-1]
            elif line.strip() and not line.startswith(("#", "//")):
                name, aligned = line.split()
                if name in seen:
                    raise ValueError(f"{path}: row {name} repeated (multi-block Stockholm)")
                seen.add(name)
                rows.append((name, aligned))
    if rf is None:
        raise ValueError(f"{path}: no #=GC RF line")
    return Path(path).name.split(".")[0], rf, rows


def summarise(path, rows_out=None):
    family, rf, rows = read_msa(path)
    model_length = sum(c != "." for c in rf)
    proteins = {}
    for name, aligned in rows:
        span = model_span(aligned, rf)
        if span is None:
            # shortcut: fail on a row with no match residue; count and skip it if real data has any
            raise ValueError(f"{path}: row {name} has no residue in a match column")
        mgyp, start, end = region(name, aligned)
        proteins.setdefault(mgyp, []).append((start, end) + span)

    s = dict.fromkeys(FIELDS, 0)
    s.update(family=family, model_length=model_length, regions=len(rows), mgyps=len(proteins))
    for mgyp, hits in proteins.items():
        copies = chain(hits, model_length)
        cls = classify(copies, model_length)
        full = [is_full(copy, model_length) for copy in copies]
        s["copies"] += len(copies)
        s["full_copies"] += sum(full)
        s["full_length_regions"] += sum(h[3] - h[2] + 1 >= FULL * model_length for h in hits)
        if cls == "single":
            s["single_full" if full[0] else "single_partial"] += 1
        else:
            s["multi_hit_mgyps"] += 1
            s[cls] += 1
            s["split_full"] += cls == "split" and full[0]
        if rows_out is not None:
            for number, (copy, copy_full) in enumerate(zip(copies, full), 1):
                for seq_from, seq_to, hmm_from, hmm_to in copy:
                    values = (
                        family,
                        mgyp,
                        seq_from,
                        seq_to,
                        hmm_from,
                        hmm_to,
                        *ends(hmm_from, hmm_to, model_length),
                        number,
                        copy_full,
                        cls,
                    )
                    rows_out.write("\t".join(map(str, values)) + "\n")
    return s


def msa_paths(arguments):
    for argument in map(Path, arguments):
        yield from sorted(argument.glob("*.sto.gz")) if argument.is_dir() else [argument]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("msa", nargs="+", help="full MSA files (.sto or .sto.gz) or folders of *.sto.gz")
    parser.add_argument("--rows", help="also write one line per region (.gz for gzip)")
    args = parser.parse_args(argv)
    rows_out = None
    if args.rows:
        rows_out = gzip.open(args.rows, "wt") if args.rows.endswith(".gz") else open(args.rows, "w")
        rows_out.write("\t".join(ROW_FIELDS) + "\n")
    totals = dict.fromkeys(FIELDS[2:], 0)
    families = 0
    print("\t".join(FIELDS))
    try:
        for path in msa_paths(args.msa):
            s = summarise(path, rows_out)
            print("\t".join(str(s[f]) for f in FIELDS), flush=True)
            families += 1
            for field in totals:
                totals[field] += s[field]
    finally:
        if rows_out is not None:
            rows_out.close()
    print(f"# total families={families} " + " ".join(f"{k}={v}" for k, v in totals.items()))


if __name__ == "__main__":
    main()
