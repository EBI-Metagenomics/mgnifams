#!/usr/bin/env python3

import argparse
import csv

def load_descriptions(desc_file):
    desc_map = {}
    with open(desc_file) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                pfam, name, description = parts[0].strip(), parts[1].strip(), parts[2].strip()
                desc_map[pfam] = {"name": name, "description": description}
            else:
                # fallback if descriptions missing, just empty strings
                pfam = parts[0]
                desc_map[pfam] = {"name": "", "description": ""}
    return desc_map

def parse_summary_block(hhr_file, desc_map):
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

                parts = line.strip().split()
                if len(parts) < 15:
                    continue
                pfam_id = parts[1]
                prob = parts[6]
                e_value = parts[7]
                length = parts[11]
                query_hmm = parts[12]
                template_hmm = parts[13] + " " + parts[14]

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
    args = parser.parse_args()

    desc_map = load_descriptions(args.descriptions)
    records = parse_summary_block(args.hhr_file, desc_map)

    with open(args.outfile, "w", newline="") as csvfile:
        fieldnames = ["id", "pfam", "name", "description", "prob", "e_value", "length", "query_hmm", "template_hmm"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            rec["id"] = args.id  # add the family chunk id to each record
            writer.writerow(rec)

    print(f"Wrote {len(records)} records to {args.outfile}")

if __name__ == "__main__":
    main()