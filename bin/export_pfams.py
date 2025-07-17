#!/usr/bin/env python3

import argparse
import csv
import re

def load_descriptions(desc_file):
    desc_map = {}
    with open(desc_file) as f:
        for line in f:
            parts = line.strip().split("\t", maxsplit=1)
            if len(parts) == 2:
                pfam, desc = parts
                desc_map[pfam] = desc
    return desc_map

def parse_summary_block(hhr_file, desc_map):
    records = []
    in_table = False

    with open(hhr_file) as f:
        for line in f:
            line = line.rstrip()

            # Detect start of No Hit table
            if line.startswith(" No Hit"):
                in_table = True
                continue

            # Blank line ends the table
            if in_table and line.strip() == "":
                in_table = False
                continue

            if in_table:
                # Parse summary table lines, example:
                # 1 PF13263.9 ; PHP_C ; PHP-associ  88.7   0.011 2.6E-06   29.4   0.0   25    2-26      4-28  (56)
                parts = re.split(r"\s+", line, maxsplit=10)
                if len(parts) < 10:
                    continue  # malformed line, skip

                pfam_full = parts[2]
                pfam = pfam_full.split(";")[0].strip()
                prob = parts[3]
                evalue = parts[4]
                cols = parts[8]
                query_hmm = parts[9]
                template_hmm = parts[10]

                desc = desc_map.get(pfam, "UNKNOWN")

                records.append({
                    "pfam": pfam,
                    "description": desc,
                    "Prob": prob,
                    "E-value": evalue,
                    "Cols": cols,
                    "Query HMM": query_hmm,
                    "Template HMM": template_hmm,
                })

    return records

def main():
    parser = argparse.ArgumentParser(description="Parse HHR summary block and output CSV with descriptions")
    parser.add_argument("--id", help="Family chunk id")
    parser.add_argument("--hhr_file", help="Input .hhr file")
    parser.add_argument("--descriptions", help="Pfam Id-Description mapping TSV file")
    parser.add_argument("--outfile", help="Output CSV file")
    args = parser.parse_args()

    desc_map = load_descriptions(args.descriptions)
    records = parse_summary_block(args.hhr_file, desc_map)

    with open(args.outfile, "w", newline="") as csvfile:
        fieldnames = ["id", "pfam", "description", "Prob", "E-value", "Cols", "Query HMM", "Template HMM"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            rec["id"] = args.id  # add the family chunk id to each record
            writer.writerow(rec)

    print(f"Wrote {len(records)} records to {args.outfile}")

if __name__ == "__main__":
    main()