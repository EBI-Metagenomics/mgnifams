#!/usr/bin/env python3

import os
import json
import argparse

def parse_tm_file(input_file, output_dir):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    assert len(lines) % 3 == 0, "Input file does not follow 3-line format per protein."

    os.makedirs(output_dir, exist_ok=True)

    for i in range(0, len(lines), 3):
        header = lines[i]
        seq = lines[i + 1]
        labels = lines[i + 2]

        if not (len(seq) == len(labels)):
            print(f"Length mismatch in entry: {header}")
            continue

        # Extract ID (e.g. >1 | GLOB → 1)
        try:
            protein_id = header.split()[0].lstrip(">")
        except:
            protein_id = f"seq_{i//3+1}"

        features = {
            "sequence": seq,
            "features": []
        }

        label_map = {
            'I': {
                "name": "Inside",
                "className": "tm-i",
                "color": "#619CFF"
            },
            'M': {
                "name": "Transmembrane",
                "className": "tm-m",
                "color": "#F8766D"
            },
            'O': {
                "name": "Outside",
                "className": "tm-o",
                "color": "#00BA38"
            }
        }

        # Group into contiguous regions
        grouped = {"I": [], "M": [], "O": []}
        current = labels[0]
        start = 0

        for j in range(1, len(labels)):
            if labels[j] != current:
                grouped[current].append({"x": start + 1, "y": j})
                start = j
                current = labels[j]
        # Append final group
        grouped[current].append({"x": start + 1, "y": len(labels)})

        # Build feature viewer JSON
        for label, rects in grouped.items():
            if rects:
                meta = label_map[label]
                features["features"].append({
                    "data": rects,
                    "name": meta["name"],
                    "className": meta["className"],
                    "color": meta["color"],
                    "type": "rect"
                })

        output_path = os.path.join(output_dir, f"{protein_id}.json")
        with open(output_path, 'w') as out:
            json.dump(features, out, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert transmembrane predictions to FeatureViewer JSON format.")
    parser.add_argument("--input_file", required=True, help="Input file with TM predictions")
    parser.add_argument("--output_dir", required=True, help="Directory to store output .json files")
    args = parser.parse_args()

    parse_tm_file(args.input_file, args.output_dir)