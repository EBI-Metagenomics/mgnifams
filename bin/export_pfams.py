#!/usr/bin/env python3

import argparse
import re
from fractions import Fraction


def load_descriptions(desc_file):
    desc_map = {}
    with open(desc_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                pfam, name, description = parts[0].strip().split(".")[0], parts[1].strip(), parts[2].strip()
                desc_map[pfam] = {"name": name, "description": description}
            else:
                # fallback if descriptions missing, just empty strings
                pfam = parts[0]
                desc_map[pfam] = {"name": "", "description": ""}
    return desc_map


def template_coverage(template_hmm):
    """`58-470 (471)` -> the fraction of the template the hit spans.

    HH-suite drops the space once the numbers fill the column: `1180-1269(1271)`.
    """
    m = re.fullmatch(r"(\d+)-(\d+)\s*\((\d+)\)", template_hmm)
    if not m:
        raise ValueError(f"unexpected Template HMM column: {template_hmm!r}")
    start, end, length = map(int, m.groups())
    return Fraction(end - start + 1, length)


def parse_summary_block(hhr_file, desc_map, min_template_coverage=0):
    records = []
    with open(hhr_file) as f:
        in_block = False
        for line in f:
            # start block after line starting with "No Hit"
            if line.startswith(" No Hit"):
                in_block = True
                continue
            if in_block:
                # empty line signals end of block
                if line.strip() == "":
                    break

                pfam_id = line[4:34].strip().split(";")[0].strip().split(".")[0]
                # Prob E-value P-value Score SS Cols Query-HMM Template-HMM: whitespace-separated, but the
                # template length may touch its span (`1180-1269(1271)`), so the rest of the line is the template
                fields = line[34:].split()
                prob, e_value, length, query_hmm = fields[0], fields[1], fields[5], fields[6]
                template_hmm = "".join(fields[7:]).replace("(", " (")

                # As FILTER_EXPORT_DOMTBL for hmmsearch: drop hits covering too little of the Pfam model
                if template_coverage(template_hmm) < min_template_coverage:
                    continue

                desc = desc_map.get(pfam_id, {"name": "", "description": ""})

                rec = {
                    "pfam": pfam_id,
                    "name": desc["name"],
                    "description": desc["description"],
                    "prob": prob,
                    "e_value": e_value,
                    "length": length,
                    "query_hmm": query_hmm,
                    "template_hmm": template_hmm,
                }
                records.append(rec)

    return records


def main():
    parser = argparse.ArgumentParser(description="Parse HHR summary block and output CSV with descriptions")
    parser.add_argument("--id", help="Family chunk id")
    parser.add_argument("--hhr_file", help="Input .hhr file")
    parser.add_argument("--descriptions", help="Pfam Id-Description mapping TSV file")
    parser.add_argument("--outfile", help="Output CSV file")
    parser.add_argument(
        "--min_template_coverage",
        type=Fraction,
        default=Fraction(0),
        help="keep a hit only if it spans at least this fraction of the template (Pfam) HMM",
    )
    args = parser.parse_args()

    desc_map = load_descriptions(args.descriptions)
    records = parse_summary_block(args.hhr_file, desc_map, args.min_template_coverage)

    with open(args.outfile, "w", encoding="utf-8") as f:
        header = ["id", "pfam", "name", "description", "prob", "e_value", "length", "query_hmm", "template_hmm"]
        f.write(",".join(header) + "\n")

        for rec in records:
            rec["id"] = args.id  # add the family chunk id to each record
            row = []
            for col in header:
                val = rec.get(col, "")
                if col == "description":
                    # Always quote description, escape any embedded quotes
                    val = '"' + val.replace('"', '""') + '"'
                row.append(val)
            f.write(",".join(row) + "\n")

    print(f"Wrote {len(records)} records to {args.outfile}")


if __name__ == "__main__":
    main()
