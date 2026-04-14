#!/usr/bin/env python3

import os
import json
import argparse
import csv

def parse_tm_file(input_file, output_dir, csv_out):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    assert len(lines) % 3 == 0, "Input file does not follow 3-line format per protein."

    os.makedirs(output_dir, exist_ok=True)

    # CSV columns for all possible labels
    csv_header = [
        "id",
        "inside_percent", "membrane_alpha_percent", "outside_percent",
        "signal_percent", "membrane_beta_percent", "periplasm_percent"
    ]
    csv_rows = [tuple(csv_header)]

    label_map = {
        'I': {
            "name": "Inside (cytosolic)",
            "className": "tm-i",
            "color": "#619CFF"
        },
        'M': {
            "name": "Transmembrane α-helix",
            "className": "tm-m",
            "color": "#F8766D"
        },
        'O': {
            "name": "Outside (extracellular)",
            "className": "tm-o",
            "color": "#00BA38"
        },
        'S': {
            "name": "Signal peptide",
            "className": "tm-s",
            "color": "#E69F00"
        },
        'B': {
            "name": "Transmembrane β-strand",
            "className": "tm-b",
            "color": "#800080"
        },
        'P': {
            "name": "Periplasm / lumen",
            "className": "tm-p",
            "color": "#00CED1"
        }
    }

    for i in range(0, len(lines), 3):
        header = lines[i]
        seq = lines[i + 1]
        labels = lines[i + 2]

        if len(seq) != len(labels):
            print(f"Length mismatch in entry: {header}")
            continue

        # Extract ID (e.g. >1 | GLOB → 1)
        try:
            protein_id = header.split()[0].lstrip(">")
        except Exception:
            protein_id = f"seq_{i//3+1}"

        features = {
            "sequence": seq,
            "features": []
        }

        total = len(labels)
        # Count percentages for each known label
        counts = {k: labels.count(k) for k in label_map.keys()}
        percents = {k: round((counts[k] / total) * 100, 2) for k in counts}

        csv_rows.append((
            protein_id,
            percents['I'], percents['M'], percents['O'],
            percents['S'], percents['B'], percents['P']
        ))

        # Group contiguous regions
        grouped = {k: [] for k in label_map.keys()}
        current = labels[0]
        start = 0

        for j in range(1, len(labels)):
            if labels[j] != current:
                if current in grouped:
                    grouped[current].append({"x": start + 1, "y": j})
                start = j
                current = labels[j]
        if current in grouped:
            grouped[current].append({"x": start + 1, "y": len(labels)})

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

    # Write CSV summary
    with open(csv_out, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(csv_rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert DeepTMHMM predictions to FeatureViewer JSON and CSV summaries.")
    parser.add_argument("--input_file", required=True, help="Input file with TM predictions")
    parser.add_argument("--output_dir", required=True, help="Directory to store output .json files")
    parser.add_argument("--csv_out", required=True, help="Path to output CSV summary file")
    args = parser.parse_args()

    parse_tm_file(args.input_file, args.output_dir, args.csv_out)
