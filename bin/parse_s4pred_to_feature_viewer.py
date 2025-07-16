#!/usr/bin/env python3

import argparse
import os
import json
import csv

def parse_s4pred_to_feature_viewer(input_file, output_file):
    conf_data, pred_data, aa_data = [], [], ""

    # Read the input file line by line
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("Conf:"):
                conf_data.extend(map(int, line.split(":")[1].strip()))
            elif line.startswith("Pred:"):
                pred_data.extend(line.split(":")[1].strip())
            elif line.startswith("  AA:"):
                aa_data += line.split(":")[1].strip()

    total = len(pred_data)
    helix_count = pred_data.count('H')
    strand_count = pred_data.count('E')
    coil_count = pred_data.count('C')
    helix_percentage = round((helix_count / total) * 100, 2)
    strand_percentage = round((strand_count / total) * 100, 2)
    coil_percentage = round((coil_count / total) * 100, 2)
    
    if not (len(conf_data) == len(pred_data) == len(aa_data)):
        raise ValueError("Mismatch in lengths of Conf, Pred, and AA sequences.")
    
    # Group Pred data into rects based on type (H - α-helices, E - β-strands, C - Coil)
    pred_features = {"H": [], "E": [], "C": []}
    start = None
    current_type = None

    for i, pred in enumerate(pred_data):
        if pred != current_type:
            if start is not None:
                pred_features[current_type].append({"x": start + 1, "y": i})
            start = i
            current_type = pred
        if i == len(pred_data) - 1:  # Handle the last region
            pred_features[current_type].append({"x": start + 1, "y": i + 1})

    # Prepare JSON for FeatureViewer
    features = {
        "sequence": aa_data,
        "features": []
    }
    
    # Add features for each Pred type (Coil, Helix, Strand)
    for key, rects in pred_features.items():
        name = {"H": "α-helices", "E": "β-strands", "C": "Coil"}.get(key, "-")
        color = {"H": "#db3089", "E": "#f2c83f", "C": "#636263"}.get(key, "#000000")
        features["features"].append({
            "data": rects,
            "name": name,
            "className": f"pred-{key.lower()}",
            "color": color,
            "type": "rect"
        })

    features["features"].append({
            "data": [{"x": i + 1, "y": conf} for i, conf in enumerate(conf_data)],
                "name": "Confidence Score",
                "className": "conf-score",
                "color": "#0F8292",
                "type": "line"
        })

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write JSON output to file
    with open(output_file, 'w') as outfile:
        json.dump(features, outfile, indent=4)

    return helix_percentage, strand_percentage, coil_percentage

def process_directory(input_dir, output_dir, csv_out):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    csv_rows = [("id","helix_percent","strand_percent","coil_percent")]

    # Iterate over all files in the input directory
    for filename in os.listdir(input_dir):
        input_file = os.path.join(input_dir, filename)

        # Only process files (skip directories)
        if os.path.isfile(input_file):
            # Construct the output file path (same basename + .json)
            base_name = os.path.splitext(filename)[0]
            output_file = os.path.join(output_dir, f"{base_name}.json")
            
            # Process the file and save the output JSON
            try:
                helix_percent, strand_percent, coil_percentage = parse_s4pred_to_feature_viewer(input_file, output_file)
                csv_rows.append((base_name, helix_percent, strand_percent, coil_percentage))
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    with open(csv_out, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(csv_rows)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse S4Pred output to JSON for FeatureViewer and collect coil % per protein.")
    parser.add_argument("--input_dir", help="Directory with S4Pred input horiz files", required=True)
    parser.add_argument("--output_dir", help="Directory for output JSON files", required=True)
    parser.add_argument("--csv_out", help="Path to output CSV summary file with coil percentages", required=True)

    args = parser.parse_args()

    process_directory(args.input_dir, args.output_dir, args.csv_out)