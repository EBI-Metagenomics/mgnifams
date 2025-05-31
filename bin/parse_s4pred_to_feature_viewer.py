#!/usr/bin/env python3

import sys
import os
import json

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

def process_directory(input_dir, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over all files in the input directory
    for filename in os.listdir(input_dir):
        input_file = os.path.join(input_dir, filename)

        # Only process files (skip directories)
        if os.path.isfile(input_file):
            # Construct the output file path (same basename + .json)
            base_name = os.path.splitext(filename)[0]
            output_file = os.path.join(output_dir, f"{base_name}.json")
            
            # Process the file and save the output JSON
            parse_s4pred_to_feature_viewer(input_file, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ./bin/parse_s4pred_to_feature_viewer.py <input_dir> <output_dir>")
        sys.exit(1)

        # Example: python ./bin/parse_s4pred_to_feature_viewer.py /path/to/s4pred_preds output/s4pred

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    process_directory(input_dir, output_dir)