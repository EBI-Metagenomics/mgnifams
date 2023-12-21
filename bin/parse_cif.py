import sys
import os
import numpy as np
from Bio.PDB import PDBParser, MMCIFIO
import subprocess

def write_header(output_file):
    header = (
        "loop_\n"
        "_ma_qa_metric.id\n"
        "_ma_qa_metric.mode\n"
        "_ma_qa_metric.name\n"
        "1 global pLDDT\n"
        "2 local  pLDDT\n"
        "#\n"
        "loop_\n"
        "_ma_qa_metric_local.label_asym_id\n"
        "_ma_qa_metric_local.label_comp_id\n"
        "_ma_qa_metric_local.label_seq_id\n"
        "_ma_qa_metric_local.metric_id\n"
        "_ma_qa_metric_local.metric_value\n"
        "_ma_qa_metric_local.model_id\n"
        "_ma_qa_metric_local.ordinal_id\n"
    )
    with open(output_file, 'w') as f:
        f.write(header)

def process_pdb(pdb_file, output_file):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line

    asym_id, comp_id, seq_id, metric_values = None, None, None, []

    with open(output_file, 'a') as f:
        for line in lines:
            if line.startswith("TER"):
                break
            parts = line.split()
            if parts[0] != "ATOM":
                continue

            new_asym_id = parts[4]
            new_comp_id = parts[3]
            new_seq_id = parts[5]

            if asym_id is not None and (new_asym_id != asym_id or new_comp_id != comp_id or new_seq_id != seq_id):
                # Write previous residue's data
                mean_metric_value = np.mean(metric_values)
                f.write(f"{asym_id} {comp_id} {seq_id} 2 {mean_metric_value:.2f} 1 1\n")
                metric_values = []

            asym_id, comp_id, seq_id = new_asym_id, new_comp_id, new_seq_id
            metric_values.append(float(parts[10]))

        # Write last residue's data
        if metric_values:
            mean_metric_value = np.mean(metric_values)
            f.write(f"{asym_id} {comp_id} {seq_id} 2 {mean_metric_value:.2f} 1 1\n")

def append_cif_model(pdb_file, output_file):
    parser = PDBParser()
    structure = parser.get_structure(os.path.splitext(pdb_file)[0], pdb_file)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save("temp2.cif")

def concat_files(file1, file2, output_file):
    command = f"cat {file1} <(tail -n +2 {file2}) > {output_file}"
    subprocess.run(command, shell=True, check=True, executable='/bin/bash')

def main():
    if len(sys.argv) != 3:
        print("Usage: python parse_to_cif.py <pdb_file> <output_file>")
        sys.exit(1)

    pdb_file    = sys.argv[1]
    output_file = sys.argv[2]
    file1       = "temp1.cif"
    file2       = "temp2.cif"

    write_header(file1)
    process_pdb(pdb_file, file1)
    append_cif_model(pdb_file, file2)
    concat_files(file1, file2, output_file)
    
    os.remove(file1)
    os.remove(file2)

if __name__ == "__main__":
    main()