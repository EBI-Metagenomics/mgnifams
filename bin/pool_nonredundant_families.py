#!/usr/bin/env python3

import argparse
import os
import shutil
import json
import csv
import gzip
import re
from concurrent.futures import ThreadPoolExecutor, as_completed


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--input_dir",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Folder of all families, including redundant ones.",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder non-redundant families.",
    )
    parser.add_argument(
        "-r",
        "--redundant",
        required=True,
        metavar="FILE",
        type=str,
        help="File with redundant family ids.",
    )
    parser.add_argument(
        "-c",
        "--similarity_csv",
        required=True,
        metavar="FILE",
        type=str,
        help="File with similarity edgelist between families.",
    )
    parser.add_argument(
        "-i",
        "--iter",
        required=True,
        metavar="INT",
        type=int,
        help="Starting family id.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        required=False,
        metavar="INT",
        type=int,
        default=None,
        help="Number of parallel threads to use.",
    )
    return parser.parse_args(args)


# mgnifam generate_families (4.x) per-chunk file names: the chunk id is the name minus its suffix
METADATA_HEADER = (
    "family_id,converged,seed_msa_size,full_msa_size,"
    "rep_protein,rep_region,rep_length,consensus_length,rep_sequence,consensus_sequence"
)
DISCARDED_HEADER = "representative,reason,value"
ARTIFACT_SUFFIXES = (".hmm.gz", ".sto.gz", ".txt")  # per-family hmm/, seed and full MSAs, rf/


def read_non_redundant_fam_ids(file_path):
    with open(file_path, "r") as file:
        return set(file.read().splitlines())


def chunk_of(filename, suffix):
    if not filename.endswith(suffix):
        raise ValueError(f"{filename}: expected a name ending in {suffix}")
    return filename[: -len(suffix)]


def split_artifact(filename):
    """`<chunk>_<n>.hmm.gz` -> (`<chunk>_<n>`, `.hmm`): the family key and the pooled (decompressed) extension."""
    for suffix in ARTIFACT_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)], suffix.removesuffix(".gz")
    raise ValueError(f"{filename}: expected one of {', '.join(ARTIFACT_SUFFIXES)}")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", buffering=1 << 20)


def create_mapping_dict():
    hmm_folder = os.path.join(arg_families_dir, "hmm")
    keys = sorted(split_artifact(name)[0] for name in os.listdir(hmm_folder))
    keys = [key for key in keys if key not in redundant_fam_ids]
    return {key: idx + arg_starting_id for idx, key in enumerate(keys)}


def pool_directory(input_dir, suffix, output_filename, splitChar, header=""):
    """Concatenates per-chunk files, rewriting each row's chunk-local family id to its pooled id. With `header`, each
    input must start with it (a changed mgnifam format would otherwise be read under the wrong column names)."""
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file = os.path.join(arg_out_dir, output_filename)

    with open(output_file, "w", buffering=1 << 20) as outfile:
        if header:
            outfile.write(header + "\n")
        for filename in sorted(os.listdir(path_to_folder)):
            chunk = chunk_of(filename, suffix)
            with open(os.path.join(path_to_folder, filename), "r", buffering=1 << 20) as infile:
                if header and infile.readline().rstrip("\n") != header:
                    raise ValueError(f"{filename}: header is not {header}")
                for line in infile:
                    value = line.strip()
                    if not value:
                        continue
                    if not splitChar:  # converged ids: one per line, families without an HMM skipped
                        fam_id = chunk + "_" + value
                        if fam_id in redundant_fam_ids or fam_id not in family_to_id:
                            continue
                        outfile.write(str(family_to_id[fam_id]) + "\n")
                        continue
                    first_element, _, remaining_elements = value.partition(splitChar)
                    fam_id = chunk + "_" + first_element
                    if fam_id in redundant_fam_ids:
                        continue
                    outfile.write(str(family_to_id[fam_id]) + splitChar + remaining_elements + "\n")


def pool_clusters_directory(input_dir, output_filename, header=""):
    """Concatenates per-chunk cluster lists; with `header`, each input starts with it and it is dropped."""
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file = os.path.join(arg_out_dir, output_filename)

    with open(output_file, "w") as outfile:
        for filename in sorted(os.listdir(path_to_folder)):
            with open(os.path.join(path_to_folder, filename), "r") as infile:
                if header and infile.readline().rstrip("\n") != header:
                    raise ValueError(f"{filename}: header is not {header}")
                content = infile.read().strip()
                if content:
                    outfile.write(content + "\n")


def translate_directory(input_dir):
    """Copies each non-redundant family's file as `<pooled id><ext>`, decompressed; HMM NAME lines get the pooled id."""
    input_folder = os.path.join(arg_families_dir, input_dir)
    output_folder = os.path.join(arg_out_dir, input_dir)
    os.makedirs(output_folder, exist_ok=True)

    for filename in sorted(os.listdir(input_folder)):
        key, extension = split_artifact(filename)
        if key in redundant_fam_ids:
            continue
        destination_path = os.path.join(output_folder, f"{family_to_id[key]}{extension}")
        with open_text(os.path.join(input_folder, filename)) as infile, open(destination_path, "w") as outfile:
            if input_dir != "hmm":
                shutil.copyfileobj(infile, outfile, 1 << 20)
                continue
            for line in infile:
                match = re.match(r"^(NAME\s+)(\S+)", line)
                if match and match.group(2) in family_to_id:
                    line = f"{match.group(1)}{family_to_id[match.group(2)]}\n"
                outfile.write(line)


def translate_edgelist(file_path, out_path):
    if os.path.getsize(file_path) == 0:
        shutil.copy(file_path, out_path)
        return

    with open(file_path, newline="") as infile, open(out_path, "w", newline="") as outfile:
        for line in infile:
            if line.startswith("#") or line.startswith("Row"):
                # Write comment or header lines as-is
                outfile.write(line)
            else:
                row = next(csv.reader([line]))
                if len(row) < 4:
                    continue  # Skip malformed lines
                row_num = row[0]
                fam1 = row[1].strip('"')
                fam2 = row[2].strip('"')
                score = row[3]

                mapped_fam1 = family_to_id.get(fam1, fam1)
                mapped_fam2 = family_to_id.get(fam2, fam2)

                # Preserve quotes in output
                outfile.write(f'{row_num},"{mapped_fam1}","{mapped_fam2}",{score}\n')


def pool_nonredundant_reps(input_folder, output_file, ids_only_file):
    path_to_folder = os.path.join(arg_families_dir, input_folder)
    output_file = os.path.join(arg_out_dir, output_file)
    ids_only_file = os.path.join(arg_out_dir, ids_only_file)

    with open(output_file, "w") as out_f, open(ids_only_file, "w") as ids_f:
        for filename in sorted(os.listdir(path_to_folder)):
            chunk_of(filename, "_reps.fasta.gz")
            with open_text(os.path.join(path_to_folder, filename)) as fasta_file:
                write_sequence = False
                for line in fasta_file:
                    if line.startswith(">"):
                        parts = line.strip().split("\t")
                        seq_id = parts[0]
                        fam_id = parts[1] if len(parts) > 1 else None

                        # Skip if description is missing or marked redundant
                        write_sequence = fam_id and fam_id not in redundant_fam_ids

                        if write_sequence:
                            # Map description to new ID
                            mapped_id = family_to_id.get(fam_id, fam_id)
                            out_f.write(f"{seq_id}\t{mapped_id}\n")
                            ids_f.write(f">{mapped_id}\n")
                    elif write_sequence:
                        out_f.write(line)
                        ids_f.write(line)


def main(args=None):
    args = parse_args(args)

    global \
        arg_families_dir, \
        arg_out_dir, \
        arg_non_redundant_fam_ids_file, \
        arg_starting_id, \
        redundant_fam_ids, \
        family_to_id

    arg_families_dir = args.input_dir
    arg_out_dir = args.output_dir
    arg_non_redundant_fam_ids_file = args.redundant
    arg_starting_id = args.iter
    arg_similarity_edgelist = args.similarity_csv
    arg_threads = args.threads

    os.makedirs(arg_out_dir, exist_ok=True)

    redundant_fam_ids = read_non_redundant_fam_ids(arg_non_redundant_fam_ids_file)
    family_to_id = create_mapping_dict()

    tasks = [
        # Same header as the update_mgnifams family_metadata.csv (both from mgnifam)
        (pool_directory, ("family_metadata", "_metadata.csv", "family_metadata.csv", ",", METADATA_HEADER)),
        (pool_directory, ("refined_families", "_families.tsv", "refined_families.tsv", "\t")),
        (pool_directory, ("converged_families", "_converged.txt", "converged_families.txt", "")),
        (pool_clusters_directory, ("successful_clusters", "successful_clusters.txt")),
        # Pooled without its header: PRESENT_DISCARDED_FAMILIES writes the MultiQC one
        (pool_clusters_directory, ("discarded_clusters", "discarded_clusters.txt", DISCARDED_HEADER)),
        (translate_directory, ("rf",)),
        (translate_directory, ("hmm",)),
        (translate_directory, ("full_msa_sto",)),
        (translate_directory, ("seed_msa_sto",)),
        (translate_edgelist, (arg_similarity_edgelist, os.path.join(arg_out_dir, "similarity_mqc.csv"))),
    ]

    with ThreadPoolExecutor(max_workers=min(len(tasks), arg_threads or os.cpu_count() or 4)) as executor:
        futures = {executor.submit(fn, *args): (fn.__name__, args) for fn, args in tasks}
        for future in as_completed(futures):
            future.result()  # re-raise any exception

    json_mapping = "family_to_id.json"
    output_file = os.path.join(arg_out_dir, json_mapping)
    with open(output_file, "w") as f:
        json.dump(family_to_id, f)

    pool_nonredundant_reps("family_reps", "family_reps.fasta", "family_ids.fasta")


if __name__ == "__main__":
    main()
