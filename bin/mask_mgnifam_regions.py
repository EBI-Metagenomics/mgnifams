#!/usr/bin/env python3
"""Mask MGnifam regions out of a (Pfam-sliced) protein FASTA, by hash buckets.

    split --kind fasta   --buckets 256 --out DIR in.fa ...        -> DIR/000.fa ... DIR/255.fa
    split --kind regions --buckets 256 --out DIR refined.tsv ...  -> DIR/000.tsv ... (protein<TAB>start<TAB>end)
    mask  --fasta DIR/007.fa --regions DIR/007.tsv --min_sequence_length 100 --out 007.masked.fa

Bucket = int(protein id) % buckets, so a protein's records and regions land in the same bucket. FASTA headers are
>{protein} (whole protein) or >{protein}_{start}_{end} (1-based inclusive slice, protein coordinates). Region rows
are family<TAB>{protein}/{start}-{end}, or family<TAB>{protein} for the whole protein (written as start 1, end 0);
1.0 names may carry trailing spaces. Every bucket file is written, empty or not.

mask cuts each record's overlapping regions out with sliceProtein, keeping protein coordinates in the headers;
a record no region overlaps is written unchanged. Masking the pieces of a Pfam-sliced FASTA gives the same pieces
as slicing the proteins on Pfam and MGnifam regions in one pass (bin/test_mask_mgnifam_regions.py).
"""

import argparse
import os
import sys

from slice_protein import sliceProtein


def as_id(text, line):
    try:
        return int(text)
    except ValueError:
        sys.exit(f"mask_mgnifam_regions: non-integer protein id {text!r} in {line!r}")


def fasta_id(header):
    """Protein id of a FASTA header line (with or without '>')."""
    return as_id(header.lstrip(">").split()[0].split("_", 1)[0], header)


def region(line):
    """(protein, start, end) of a refined_families row; end 0 = to the protein's end."""
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 2:
        sys.exit(f"mask_mgnifam_regions: region row without a member column: {line!r}")
    member = cols[1].strip()
    protein, _, span = member.partition("/")
    if not span:
        return as_id(protein, line), 1, 0
    start, _, end = span.partition("-")
    return as_id(protein, line), as_id(start, line), as_id(end, line)


def split(kind, buckets, out, inputs):
    os.makedirs(out, exist_ok=True)
    ext = "fa" if kind == "fasta" else "tsv"
    handles = [open(os.path.join(out, f"{b:03d}.{ext}"), "w", buffering=1 << 20) for b in range(buckets)]
    try:
        for path in inputs:
            with open(path, buffering=1 << 20) as fh:
                if kind == "fasta":
                    handle = None
                    for line in fh:
                        if line.startswith(">"):
                            handle = handles[fasta_id(line) % buckets]
                        elif handle is None:
                            sys.exit(f"mask_mgnifam_regions: {path} does not start with a FASTA header")
                        handle.write(line)
                else:
                    for line in fh:
                        if line.strip():
                            protein, start, end = region(line)
                            handles[protein % buckets].write(f"{protein}\t{start}\t{end}\n")
    finally:
        for h in handles:
            h.close()


def read_fasta(path):
    header, chunks = None, []
    with open(path, buffering=1 << 20) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def mask_record(header, sequence, regions, min_sequence_length):
    """FASTA records left of one record after removing the regions (protein coordinates; end 0 = protein end)."""
    parts = header.split("_")
    if len(parts) not in (1, 3):
        sys.exit(f"mask_mgnifam_regions: header >{header} is neither >protein nor >protein_start_end")
    protein = parts[0]
    start, end = (int(parts[1]), int(parts[2])) if len(parts) == 3 else (1, len(sequence))
    if end - start + 1 != len(sequence):
        sys.exit(f"mask_mgnifam_regions: >{header} spans {end - start + 1} residues, its sequence {len(sequence)}")
    local = []
    for s, e in regions:
        e = end if e == 0 else e
        if s <= end and e >= start:
            local.append(["m", max(s, start) - start + 1, min(e, end) - start + 1])
    if not local:
        return [f">{header}\n{sequence}\n"]
    return sliceProtein(protein, sequence, {"m": local}, min_sequence_length, start=start)


def mask(fasta, regions_path, min_sequence_length, out):
    regions = {}
    with open(regions_path, buffering=1 << 20) as fh:
        for line in fh:
            protein, start, end = line.split("\t")
            regions.setdefault(int(protein), []).append((int(start), int(end)))
    with open(out, "w", buffering=1 << 20) as fh:
        for header, sequence in read_fasta(fasta):
            fh.writelines(mask_record(header, sequence, regions.get(fasta_id(header), ()), min_sequence_length))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("split")
    p.add_argument("--kind", choices=("fasta", "regions"), required=True)
    p.add_argument("--buckets", type=int, default=256)
    p.add_argument("--out", required=True)
    p.add_argument("inputs", nargs="+")
    m = sub.add_parser("mask")
    m.add_argument("--fasta", required=True)
    m.add_argument("--regions", required=True)
    m.add_argument("--min_sequence_length", type=int, required=True)
    m.add_argument("--out", required=True)
    a = parser.parse_args(argv)
    if a.command == "split":
        split(a.kind, a.buckets, a.out, a.inputs)
    else:
        mask(a.fasta, a.regions, a.min_sequence_length, a.out)


if __name__ == "__main__":
    main()
